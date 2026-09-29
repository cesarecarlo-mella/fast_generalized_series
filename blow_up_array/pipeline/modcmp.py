#!/usr/bin/env python3
"""
modcmp.py -- compare two Mathematica list files (e.g. phys_*.m / ber_*.m from
two different pipelines) component by component, WITHOUT expanding anything:
every component is evaluated modulo a 61-bit prime p with each symbol
(x, y, Log[x], Pi, Zeta[3], sx, LogY, I, ...) replaced by a fixed pseudo-random
residue.  Two expressions that differ as polynomials/Laurent polynomials agree
only with probability ~ degree/p (~1e-16): a mismatch is a real difference,
a match is equality with overwhelming probability (Schwartz-Zippel).

    python3 modcmp.py fileA.m fileB.m            -> prints IDENTICAL / DIFFERENT
Pure Python 3.8+, no packages.  Streams tokens, so memory stays ~file size.
"""
import hashlib
import re
import sys
from fractions import Fraction


def _is_prime(n):
    if n < 2:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % q == 0:
            return n == q
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


P = (1 << 61) - 1
while not (P % 4 == 1 and _is_prime(P)):   # p = 1 mod 4 so that I = sqrt(-1) exists
    P -= 2
for a in range(2, 1000):
    c = pow(a, (P - 1) // 4, P)
    if c * c % P == P - 1:
        SQRT_M1 = c
        break

_root = {}


def root_of(name):
    """r_s: the value of Sqrt[s]; the symbol itself is r_s^2."""
    r = _root.get(name)
    if r is None:
        h = hashlib.sha256(("modcmp-v1:" + name).encode()).digest()
        r = int.from_bytes(h[:16], 'big') % P or 7
        _root[name] = r
    return r


def sym_value(name):
    if name == 'I':
        return SQRT_M1
    r = root_of(name)
    return r * r % P


_TOK = re.compile(r'\s*(?:(\d+)|([A-Za-z][A-Za-z0-9$]*)((?:\[[^\[\]]*\])*)|(->|[-+*/^(){},]))')


POS = [0]      # end offset of the last token handed out (for resuming)


def tokens(s, pos=0):
    for m in _TOK.finditer(s, pos):
        num, sym, br, op = m.groups()
        POS[0] = m.end()
        if num is not None:
            yield ('n', num)
        elif sym is not None:
            yield ('s', sym + (re.sub(r'\s*,\s*', ', ', re.sub(r'\s+', '', br)) if br else ''))
        elif op is not None:
            yield ('o', op)
    yield ('e', None)


class Ev:
    """Pratt evaluator mod P.  Values are ints mod P; an atom that is a bare
    symbol also carries its name so fractional powers (Sqrt, x^(3/2)) work."""

    def __init__(self, s, pos=0):
        self.it = tokens(s, pos)
        self.tok = next(self.it)

    def adv(self):
        self.tok = next(self.it)

    def expect(self, op):
        if self.tok != ('o', op):
            raise ValueError("expected %s, got %r" % (op, self.tok))
        self.adv()

    def top(self):
        """a top-level list -> list of residues (nested lists flattened in order)."""
        self.expect('{')
        out = []
        if self.tok != ('o', '}'):
            out.append(self.item())
            while self.tok == ('o', ','):
                self.adv(); out.append(self.item())
        self.expect('}')
        return out

    def item(self):
        if self.tok == ('o', '{'):          # nested list: fold into one residue deterministically
            vals = self.top()
            acc = 0
            for k, v in enumerate(vals):
                acc = (acc * 1000003 + v + k) % P
            return acc
        return self.expr()

    def expr(self):
        neg = False
        if self.tok == ('o', '-'):
            self.adv(); neg = True
        elif self.tok == ('o', '+'):
            self.adv()
        v = self.term()
        if neg:
            v = -v % P
        while True:
            t = self.tok
            if t == ('o', '+'):
                self.adv(); v = (v + self.term()) % P
            elif t == ('o', '-'):
                self.adv(); v = (v - self.term()) % P
            else:
                return v

    def term(self):
        v = self.unary()
        while True:
            t = self.tok
            if t == ('o', '*'):
                self.adv(); v = v * self.unary() % P
            elif t == ('o', '/'):
                self.adv(); v = v * pow(self.unary(), P - 2, P) % P
            else:
                return v

    def unary(self):
        if self.tok == ('o', '-'):
            self.adv()
            return -self.unary() % P
        return self.power()

    def power(self):
        v, name = self.atom()
        if self.tok == ('o', '^'):
            self.adv()
            e = self.exponent()
            if e.denominator == 1:
                n = int(e)
                return pow(v, n, P) if n >= 0 else pow(pow(v, P - 2, P), -n, P)
            if e.denominator == 2 and name is not None and name != 'I':
                k = int(e * 2)
                r = root_of(name)
                return pow(r, k, P) if k >= 0 else pow(pow(r, P - 2, P), -k, P)
            raise ValueError("unsupported power %s of %r" % (e, name))
        return v

    def exponent(self):
        """exponents are small exact rationals: 3, (-1), (3/2), -2 ..."""
        sign = 1
        if self.tok == ('o', '-'):
            self.adv(); sign = -1
        if self.tok[0] == 'n':
            e = Fraction(int(self.tok[1])); self.adv(); return sign * e
        self.expect('(')
        num = []
        depth = 1
        while depth:
            t = self.tok
            if t == ('o', '('):
                depth += 1
            elif t == ('o', ')'):
                depth -= 1
                if depth == 0:
                    break
            num.append(t[1]); self.adv()
        self.adv()
        return sign * Fraction(eval("".join(num), {"__builtins__": {}}, {})).limit_denominator(10 ** 6) \
            if not any(c.isalpha() for c in "".join(num)) else self._bad("".join(num))

    def _bad(self, s):
        raise ValueError("non-numeric exponent %s" % s)

    def atom(self):
        kind, val = self.tok
        if kind == 'n':
            self.adv(); return int(val) % P, None
        if kind == 's':
            self.adv()
            if val.startswith('Sqrt['):
                inner = val[5:-1]
                if re.fullmatch(r'[A-Za-z]\w*', inner):
                    return root_of(inner), None
                raise ValueError("Sqrt of non-symbol: %s" % val)
            return sym_value(val), (val if '[' not in val else None)
        if val == '(':
            self.adv()
            v = self.expr()
            self.expect(')')
            return v, None
        raise ValueError("unexpected token %r" % ((kind, val),))


def evaluate_file(path):
    with open(path) as f:
        s = f.read()
    s = re.sub(r'\(\*.*?\*\)', ' ', s, flags=re.S)
    s = re.sub(r'\\\r?\n\s*', '', s)
    return Ev(s).top()


def eval_resumable(path, state, budget):
    """evaluate the top-level list in `path` for at most `budget` seconds,
    continuing from / saving to the json file `state`.  Returns True when done."""
    import json, os, time
    t0 = time.time()
    st = json.load(open(state)) if os.path.isfile(state) else {"pos": None, "vals": [], "done": False}
    if st["done"]:
        return True
    s = _prep(open(path).read())
    if st["pos"] is None:                       # skip the opening '{'
        st["pos"] = s.index('{') + 1
    ev = Ev(s, st["pos"])
    if ev.tok == ('o', '}'):
        st["done"] = True
    while not st["done"] and time.time() - t0 < budget:
        st["vals"].append(ev.item())
        t = ev.tok                              # ',' or '}'
        st["pos"] = POS[0]                      # just after that separator
        if t == ('o', '}'):
            st["done"] = True
        elif t != ('o', ','):
            raise ValueError("bad separator %r" % (t,))
        else:
            ev.adv()                            # tokenizer now points at the next item
            ev = Ev(s, st["pos"])               # (re-anchor so POS/state stay consistent)
    json.dump(st, open(state + ".tmp", "w")); os.replace(state + ".tmp", state)
    return st["done"]


def _prep(s):
    s = re.sub(r'\(\*.*?\*\)', ' ', s, flags=re.S)
    return re.sub(r'\\\r?\n\s*', '', s)


if __name__ == '__main__' and len(sys.argv) > 1 and sys.argv[1] == '--resume':
    # modcmp.py --resume FILE STATE BUDGET_SECONDS  -> prints DONE or PARTIAL n
    done = eval_resumable(sys.argv[2], sys.argv[3], float(sys.argv[4]))
    import json
    print(("DONE " if done else "PARTIAL ") + str(len(json.load(open(sys.argv[3]))["vals"])), sys.argv[2])
    sys.exit(0)

if __name__ == '__main__':
    a, b = sys.argv[1], sys.argv[2]
    va, vb = evaluate_file(a), evaluate_file(b)
    if len(va) != len(vb):
        print("DIFFERENT  %s: %d vs %d components" % (b.split('/')[-1], len(va), len(vb))); sys.exit(1)
    bad = [i + 1 for i in range(len(va)) if va[i] != vb[i]]
    name = a.split('/')[-1]
    if bad:
        print("DIFFERENT  %-32s %d/%d components differ, e.g. %s" % (name, len(bad), len(va), bad[:10]))
        sys.exit(1)
    print("IDENTICAL  %-32s %d components" % (name, len(va)))
