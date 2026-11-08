"""Reproduce the explicitly partial native-Windows Jacobian 0.23.0 runtime.

The official wheel is unchanged. PPLpy is unavailable on this platform; its
absence is reported, never patched out of Jacobian's distribution metadata.
"""
from pathlib import Path
from datetime import datetime,timezone
import argparse,hashlib,json,os,platform,shutil,subprocess,sys,urllib.request,venv

VERSION='0.23.0'
WHEEL_SHA256='72fe5a9210db17e18e756a25cb94473b2c188bba591ad8ab8d7dfff9531e05f8'
HERE=Path(__file__).resolve().parent


def need(ok,label):
    if not ok:raise RuntimeError(label)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--register-mcp',action='store_true',help='register this pinned runtime with the Codex CLI')
    parser.add_argument('--dry-run',action='store_true')
    args=parser.parse_args()
    need(sys.platform=='win32' and (3,12)<=sys.version_info[:2]<(3,14),'requires 64-bit Windows CPython 3.12 or 3.13')
    need(platform.architecture()[0]=='64bit','requires 64-bit Python')
    runtime=Path.home()/'.codex/tools/jacobian'/VERSION;py=runtime/'venv/Scripts/python.exe'
    requirements=[q for q in (HERE/'requirements-windows.lock').read_text().splitlines()
                  if q.strip() and not q.lower().startswith('jacobian==')]
    need(not any(q.lower().startswith('pplpy') for q in requirements),'explicit native-Windows dependency scope')
    if args.dry_run:
        print(json.dumps({'version':VERSION,'runtime':str(runtime),'register_mcp':args.register_mcp,
            'unavailable':['pplpy','Singular','QEPCAD'],'writes_performed':False},indent=2));return
    runtime.mkdir(parents=True,exist_ok=True)
    if not py.exists():venv.EnvBuilder(with_pip=True).create(runtime/'venv')
    env=os.environ.copy()
    for key in ('http','https'):
        proxy=urllib.request.getproxies().get(key)
        if proxy:env[key.upper()+'_PROXY']=proxy
    # urllib uses the platform's configured proxy. Never disable TLS checks.
    release=json.load(urllib.request.urlopen('https://pypi.org/pypi/jacobian/'+VERSION+'/json',timeout=30))
    item=next(q for q in release['urls'] if q['filename']=='jacobian-'+VERSION+'-py3-none-any.whl')
    need(item['digests']['sha256']==WHEEL_SHA256,'fixed official release digest')
    wheel=runtime/item['filename']
    raw=urllib.request.urlopen(item['url'],timeout=60).read()
    need(hashlib.sha256(raw).hexdigest()==WHEEL_SHA256,'downloaded official wheel digest')
    wheel.write_bytes(raw)
    subprocess.run([str(py),'-m','pip','install','--only-binary=:all:',*requirements],env=env,check=True)
    subprocess.run([str(py),'-m','pip','install','--no-deps',str(wheel)],env=env,check=True)
    check=subprocess.run([str(py),'-m','pip','check'],capture_output=True,text=True)
    issues=[q for q in check.stdout.splitlines() if q.strip()]
    need(not check.returncode or (bool(issues) and all('pplpy' in q.lower() for q in issues)),
         'unexpected dependency conflict: '+check.stdout+check.stderr)
    subprocess.run([str(py),'-m','jacobian.mcp.stdio_cli','--version'],check=True)
    registered=False
    if args.register_mcp:
        cli=shutil.which('codex.cmd') or shutil.which('codex');need(cli is not None,'Codex CLI missing')
        result=subprocess.run([cli,'mcp','list','--json'],capture_output=True,text=True,encoding='utf-8',check=True)
        old=next((s for s in json.loads(result.stdout) if s['name']=='jacobian'),None)
        if old:
            transport=old['transport']
            need(Path(transport['command']).resolve()==py.resolve() and transport['args']==['-m','jacobian.mcp.stdio_cli'],
                 'existing foreign Jacobian registration preserved; inspect it before replacement')
            need(old['enabled'],'existing matching registration is disabled; inspect before changing')
        else:
            config_root=Path(os.environ.get('CODEX_HOME',str(Path.home()/'.codex')))
            config=config_root/'config.toml'
            if config.exists():
                stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
                shutil.copy2(config,runtime/('config-before-jacobian-'+stamp+'.toml'))
            subprocess.run([cli,'mcp','add','jacobian','--env','PYTHONUTF8=1','--',str(py),'-m','jacobian.mcp.stdio_cli'],check=True)
        registered=True
    print(json.dumps({'installed_version':VERSION,'official_wheel_sha256':WHEEL_SHA256,
        'MCP_registered':registered,'pplpy_enabled':False,'dependency_check_exceptions':issues},indent=2))


if __name__=='__main__':main()
