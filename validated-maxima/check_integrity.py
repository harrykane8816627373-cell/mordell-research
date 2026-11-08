"""Reject altered Newton/Hessian claims and altered fixed bytes; -O safe."""
from pathlib import Path
from copy import deepcopy
import json, shutil, subprocess, sys, tempfile
from verify_all import verify_row

HERE = Path(__file__).resolve().parent


def main():
    common = json.loads((HERE/'certificates.json').read_text(encoding='utf-8'))
    stable = next(row for row in common['rows'] if row['n'] == 10)
    saddle = next(row for row in common['rows'] if row['n'] == 15)
    cases = []
    row = deepcopy(stable)
    row['inverse_numerators'] = [[0 for q in line] for line in row['inverse_numerators']]
    cases.append(('zero preconditioner',row))
    row = deepcopy(stable)
    row['whitener_numerators'] = [[0 for q in line] for line in row['whitener_numerators']]
    cases.append(('singular alleged Hessian congruence',row))
    row = deepcopy(stable)
    row['center_numerators'][0] += 1 << (common['center_exponent']-3)
    cases.append(('shift midpoint outside root box',row))
    row = deepcopy(saddle)
    row['positive_vector_numerators'] = [0 for q in row['positive_vector_numerators']]
    cases.append(('zero alleged instability direction',row))
    rejected = []
    for label,row in cases:
        try:
            verify_row(row,common)
        except ArithmeticError:
            rejected.append(label)
        else:
            raise ArithmeticError('false certificate accepted: '+label)
    with tempfile.TemporaryDirectory(prefix='mordell-local-tamper-') as td:
        stage = Path(td)
        for path in HERE.iterdir():
            if path.is_file():
                shutil.copyfile(path,stage/path.name)
        source = stage/'certificates.json'
        source.write_bytes(source.read_bytes()+b'\n')
        proc = subprocess.run([sys.executable,'-B','-S','-O',str(stage/'verify_all.py')],
            cwd=stage,capture_output=True,text=True,encoding='utf-8',timeout=30)
        if proc.returncode == 0 or 'fixed input certificates.json' not in proc.stderr:
            raise ArithmeticError('fixed-byte tamper was not rejected')
        rejected.append('one altered byte in fixed certificate file')
    result = {'passed':True,'optimization_safe':True,'false_certificates_rejected':rejected}
    (HERE/'integrity_report.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
