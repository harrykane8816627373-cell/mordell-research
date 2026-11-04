"""Bind the already accepted 19dbbd3 proof and its unchanged arithmetic inputs."""
from pathlib import Path,PurePosixPath
import gzip,hashlib,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SHA256='6e7a0b073e8ea28f0a2aceefe2b1a466cbf118c3e2475f46f57ab0417d96b081'
META_SHA256='9e3a8762b8c0bab9fecc9efebfc07312963460f2ae4bad5e774d7a57336a465c'
IDENTITY_KEY='nov2_prefix_sha256'


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def load_current():
    raw=(HERE/'certificates/force_checkpoint_19dbbd3.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest()==META_SHA256,'fixed accepted 19dbbd3 metadata')
    meta=json.loads(raw)
    raw=(HERE/'certificates/force_checkpoint_19dbbd3.json.gz').read_bytes()
    need(meta['certificate_sha256']==SHA256 and hashlib.sha256(raw).hexdigest()==SHA256,'fixed accepted 19dbbd3 certificate')
    plain=gzip.decompress(raw)
    need(hashlib.sha256(plain).hexdigest()==meta['decompressed_sha256'],'accepted 19dbbd3 decompressed identity')
    base=json.loads(plain)
    raw=(HERE/'certificates/force_checkpoint_19dbbd3_report.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest()==meta['report_sha256'],'fixed accepted 19dbbd3 report')
    report=json.loads(raw)
    need(report['passed'] is True and report['force_certificate_sha256']==SHA256
        and report['processed_nodes']==meta['processed_nodes']==94399
        and report['open_leaf_boxes']==meta['pending_leaves']==14
        and report['full_eight_point_exclusion_proved'] is False,'accepted incomplete 19dbbd3 prefix')
    names=set()
    for row in meta['arithmetic_inputs']:
        name=row['path'];p=PurePosixPath(name)
        need(not p.is_absolute() and '..' not in p.parts and ':' not in name and '\\' not in name
            and p.as_posix()==name and name not in names,'safe unique Nov 2 arithmetic input')
        names.add(name);raw=(ROOT/name).read_bytes()
        need(len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256'],
            'unchanged accepted 19dbbd3 arithmetic '+name)
    need(len(names)==meta['arithmetic_input_count']==70,'all Nov 2 arithmetic inputs checked')
    return base,report,meta


def descendant(data,base):
    need(data.get(IDENTITY_KEY)==SHA256,'explicit accepted 19dbbd3 prefix identity')
    for key in ['parent_sha256','attachment_sha256','precision_bits','scope','kernel_version',
            'checkpoint_sha256','accepted_prefix_sha256','roots']:
        need(data[key]==base[key],'unchanged Nov 2 prefix model '+key)
    need(len(data['nodes'])>=len(base['nodes']),'complete Nov 2 prefix retained')
    for i,node in enumerate(base['nodes']):
        if node['kind']!='unknown':need(data['nodes'][i]==node,'every accepted 19dbbd3 step is unchanged')
