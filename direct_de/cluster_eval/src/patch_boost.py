"""copy boost's odeint bulirsch_stoer.hpp into BUILD/boostpatch/... and make it compile with QD's dd_real:
explicit double casts of the integer step-sequence values (size_t -> dd_real is ambiguous).
The algorithm is unchanged.   usage: patch_boost.py BUILD [BOOST_INC]"""
import os, re, subprocess, sys

build = sys.argv[1]
inc = sys.argv[2] if len(sys.argv) > 2 and sys.argv[2] else None
rel = 'boost/numeric/odeint/stepper/bulirsch_stoer.hpp'
cands = [inc] if inc else []
if not inc:                              # compiler's default include path
    out = subprocess.run(['sh', '-c', '${CXX:-g++} -E -x c++ - -v < /dev/null'], capture_output=True, text=True).stderr
    on = False
    for line in out.splitlines():
        if line.startswith('#include <...>'): on = True; continue
        if line.startswith('End of search list'): break
        if on: cands.append(line.strip())
    cands += ['/usr/include', '/usr/local/include']
src = next((os.path.join(c, rel) for c in cands if c and os.path.exists(os.path.join(c, rel))), None)
if not src:
    sys.exit('boost odeint not found (set BOOST_INC to the directory that contains boost/)')
s = open(src).read()
n0 = len(s)
s = re.sub(r'static_cast< ?value_type ?>\( ?(m_interval_sequence\[[a-z_+0-9]+\]|m_cost\[[a-z_+0-9]+\]) ?\)',
           r'static_cast<value_type>( static_cast<double>(\1) )', s)
s = re.sub(r'const value_type d = (m_interval_sequence\[m_current_k_opt\] \* m_interval_sequence\[m_current_k_opt\+1\] /\s*\(m_interval_sequence\[0\]\*m_interval_sequence\[0\]\));',
           r'const value_type d = static_cast<double>( \1 );', s)
s = re.sub(r'const value_type d = (m_interval_sequence\[m_current_k_opt\] / m_interval_sequence\[0\]);',
           r'const value_type d = static_cast<double>( \1 );', s)
dst = os.path.join(build, 'boostpatch', rel)
os.makedirs(os.path.dirname(dst), exist_ok=True)
open(dst, 'w').write(s)
print('patched', src, '->', dst, '(%d casts)' % s.count('static_cast<double>'))
