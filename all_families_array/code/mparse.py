"""
mparse.py -- minimal, exact parser for the Mathematica InputForm files used by
the blow-up pipeline (atilde.m, dlog_d*_c*.m, sol0.m, bdy_*.m, exact_*.m).

Expressions are parsed into sparse multivariate "polynomials" over the
rationals:  Poly = dict  key -> gmpy2.mpq,  key = sorted tuple of
(atom, exp2) with exp2 = 2*exponent (so half-integer powers such as
Sqrt[t] or t^(3/2) are integers here).  The imaginary unit I is an atom
whose exponent is reduced with I^2 = -1, so Gaussian rationals are exact.

Supported grammar (everything that appears in the pipeline files):
    numbers, symbols, Head[...] atoms (Log[t], Zeta[3], l[5], Sqrt[t]),
    + - * / ^ ( ) { } , ->
Division is only allowed by a single-term Poly (a number or a monomial),
powers of sums only with non-negative integer exponents -- both are all
the files ever need.
"""
import re
import gmpy2
from gmpy2 import mpq

# a symbol may carry several bracket groups, e.g. the unresolved boundary
# constants  cnst[6, 13][5]  -- kept as ONE opaque atom
_TOK = re.compile(r'\s*(?:(\d+)|([A-Za-z][A-Za-z0-9$]*)((?:\[[^\[\]]*\])*)|(->|[-+*/^(){},]))')
I_ATOM = 'I'
ONE = ()


def _mulkey(a, b):
    """multiply two monomial keys; returns (key, sign)."""
    if not a:
        return b, 1
    if not b:
        return a, 1
    d = dict(a)
    for at, e in b:
        d[at] = d.get(at, 0) + e
    sign = 1
    ei = d.get(I_ATOM)
    if ei is not None:
        # I exponent in units of 2 (exp2); I^2 -> -1
        k = ei // 2
        r = k % 4
        if r == 2 or r == 3:
            sign = -1
        if r % 2 == 1:
            d[I_ATOM] = 2
        else:
            del d[I_ATOM]
    return tuple(sorted((at, e) for at, e in d.items() if e != 0)), sign


def padd(a, b, sb=1):
    if len(a) < len(b) and sb == 1:
        a, b = b, a
    r = dict(a)
    for k, v in b.items():
        nv = r.get(k, 0) + (v if sb == 1 else -v)
        if nv:
            r[k] = nv
        elif k in r:
            del r[k]
    return r


def pmul(a, b):
    if len(a) == 1 and ONE in a:
        c = a[ONE]
        return {k: v * c for k, v in b.items()} if c else {}
    if len(b) == 1 and ONE in b:
        c = b[ONE]
        return {k: v * c for k, v in a.items()} if c else {}
    r = {}
    for ka, va in a.items():
        for kb, vb in b.items():
            k, s = _mulkey(ka, kb)
            v = va * vb if s == 1 else -va * vb
            nv = r.get(k, 0) + v
            if nv:
                r[k] = nv
            elif k in r:
                del r[k]
    return r


def pconst(p):
    """value of a constant Poly (or None)."""
    if not p:
        return mpq(0)
    if len(p) == 1 and ONE in p:
        return p[ONE]
    return None


def pinv_monomial(p):
    if len(p) != 1:
        raise ValueError("division by a non-monomial: %r" % (p,))
    (k, v), = p.items()
    kk = tuple((at, -e) for at, e in k)
    # I^-1 = -I
    sign = 1
    kk2 = []
    for at, e in kk:
        if at == I_ATOM:
            kk2.append((at, 2)); sign = -1
        else:
            kk2.append((at, e))
    return {tuple(kk2): sign / v}


def ppow(base, ex):
    c = pconst(ex)
    if c is None:
        raise ValueError("non-constant exponent")
    if len(base) == 1:
        (k, v), = base.items()
        if k == ONE:
            if c.denominator != 1:
                raise ValueError("fractional power of a number")
            return {ONE: v ** int(c)}
        if v == 1 and not any(at == I_ATOM for at, _ in k):
            e2 = [(at, e * c) for at, e in k]
            for at, e in e2:
                if e.denominator != 1:
                    raise ValueError("exponent not representable in half-units")
            return {tuple((at, int(e)) for at, e in e2): mpq(1)}
    if c.denominator != 1 or c < 0:
        raise ValueError("power of a sum with exponent %s" % c)
    r = {ONE: mpq(1)}
    for _ in range(int(c)):
        r = pmul(r, base)
    return r


class Parser:
    def __init__(self, s):
        s = re.sub(r'\(\*.*?\*\)', ' ', s, flags=re.S)      # comments
        s = re.sub(r'\\\r?\n\s*', '', s)                      # line continuation
        self.toks = []
        pos = 0
        n = len(s)
        app = self.toks.append
        for m in _TOK.finditer(s):
            if m.start() != pos and s[pos:m.start()].strip():
                raise ValueError("tokenizer stuck at %r" % s[pos:pos + 40])
            pos = m.end()
            num, sym, br, op = m.groups()
            if num is not None:
                app(('n', num))
            elif sym is not None:
                app(('s', sym + re.sub(r'\s*,\s*', ', ', re.sub(r'\s+', '', br or ''))))
            else:
                app(('o', op))
        if s[pos:].strip():
            raise ValueError("trailing garbage %r" % s[pos:pos + 40])
        self.i = 0

    def peek(self):
        return self.toks[self.i] if self.i < len(self.toks) else ('e', None)

    def take(self):
        t = self.toks[self.i]; self.i += 1; return t

    def expect(self, op):
        t = self.take()
        if t != ('o', op):
            raise ValueError("expected %s got %r" % (op, t))

    # ---- grammar -------------------------------------------------------
    def parse(self):
        v = self.value()
        if self.i != len(self.toks):
            raise ValueError("unparsed tail at token %d: %r" % (self.i, self.toks[self.i:self.i + 5]))
        return v

    def value(self):
        """list, rule, or expression."""
        if self.peek() == ('o', '{'):
            self.take()
            out = []
            if self.peek() != ('o', '}'):
                out.append(self.value())
                while self.peek() == ('o', ','):
                    self.take(); out.append(self.value())
            self.expect('}')
            return out
        e = self.expr()
        if self.peek() == ('o', '->'):
            self.take()
            return ('rule', e, self.value())
        return e

    def expr(self):
        if self.peek() == ('o', '-'):
            self.take(); r = self.term(); r = {k: -v for k, v in r.items()}
        else:
            if self.peek() == ('o', '+'):
                self.take()
            r = self.term()
        while True:
            t = self.peek()
            if t == ('o', '+'):
                self.take(); r = padd(r, self.term())
            elif t == ('o', '-'):
                self.take(); r = padd(r, self.term(), -1)
            else:
                return r

    def term(self):
        r = self.unary()
        while True:
            t = self.peek()
            if t == ('o', '*'):
                self.take(); r = pmul(r, self.unary())
            elif t == ('o', '/'):
                self.take(); r = pmul(r, pinv_monomial(self.unary()))
            else:
                return r

    def unary(self):
        if self.peek() == ('o', '-'):
            self.take()
            r = self.unary()
            return {k: -v for k, v in r.items()}
        return self.power()

    def power(self):
        b = self.atom()
        if self.peek() == ('o', '^'):
            self.take()
            e = self.unary()
            return ppow(b, e)
        return b

    def atom(self):
        kind, v = self.take()
        if kind == 'n':
            iv = int(v)
            return {ONE: mpq(iv)} if iv else {}
        if kind == 's':
            if v.startswith('Sqrt['):
                inner = v[5:-1].strip()
                if not re.fullmatch(r'[A-Za-z]\w*', inner):
                    raise ValueError("Sqrt of non-symbol: %s" % v)
                return {((inner, 1),): mpq(1)}
            if v == 'I':
                return {((I_ATOM, 2),): mpq(1)}
            return {((v, 2),): mpq(1)}
        if v == '(':
            e = self.expr()
            self.expect(')')
            return e
        raise ValueError("unexpected token %r" % ((kind, v),))


def parse_string(s):
    return Parser(s).parse()


def parse_file(path):
    with open(path) as f:
        return parse_string(f.read())


def split_top_list(s):
    """split a top-level Mathematica list {a, b, ...} into element strings
    (lets huge files be parsed element by element / in parallel)."""
    s = re.sub(r'\(\*.*?\*\)', ' ', s, flags=re.S)
    s = s.strip()
    assert s[0] == '{' and s[-1] == '}', "not a list"
    out, depth, start = [], 0, 1
    for idx in range(1, len(s) - 1):
        c = s[idx]
        if c in '({[':
            depth += 1
        elif c in ')}]':
            depth -= 1
        elif c == ',' and depth == 0:
            out.append(s[start:idx]); start = idx + 1
    out.append(s[start:len(s) - 1])
    return out


# ---------------------------------------------------------------- writer
def fmt_rat(q):
    return str(q.numerator) if q.denominator == 1 else "%d/%d" % (q.numerator, q.denominator)


def fmt_poly(p, order=None):
    """InputForm-ish string for a Poly (used to write exact_*.m files)."""
    if not p:
        return "0"
    terms = []
    for k, v in (sorted(p.items(), key=order) if order else p.items()):
        fac = []
        for at, e in k:
            if at == I_ATOM:
                fac.append("I")
            elif e == 2:
                fac.append(at)
            elif e % 2 == 0:
                fac.append("%s^%d" % (at, e // 2) if e > 0 else "%s^(%d)" % (at, e // 2))
            else:
                fac.append("%s^(%d/2)" % (at, e))
        mono = "*".join(fac)
        if not mono:
            terms.append("(" + fmt_rat(v) + ")")
        else:
            terms.append("(" + fmt_rat(v) + ")*" + mono)
    return " + ".join(terms)
