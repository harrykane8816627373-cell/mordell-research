"""Validate an explicitly reused, previously replayed force-proof checkpoint.

This checks identity and immutable arithmetic inputs, not the old arithmetic
proof. The full checker remains available without checkpoint reuse.
"""
from pathlib import Path,PurePosixPath
import gzip,hashlib,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SHA256='2f2af13b0164ded1fc5b199874c76f790b0e3c713fa4b432b36d86525d22b30b'
META_SHA256='213ccb3f16f35fd8067255ffb8d4b8e659c47e7a15e66d9b75c00f4211087d8c'


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def load_checkpoint():
    raw=(HERE/'certificates/force_checkpoint_768598a.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest()==META_SHA256,'fixed accepted checkpoint metadata')
    meta=json.loads(raw)
    raw=(HERE/'certificates/force_checkpoint_768598a.json.gz').read_bytes()
    need(meta['certificate_sha256']==SHA256 and hashlib.sha256(raw).hexdigest()==SHA256,
        'fixed previously replayed checkpoint bytes')
    plain=gzip.decompress(raw)
    need(hashlib.sha256(plain).hexdigest()==meta['decompressed_sha256'],'checkpoint decompressed identity')
    base=json.loads(plain)
    report_raw=(HERE/'certificates/force_checkpoint_768598a_report.json').read_bytes()
    need(hashlib.sha256(report_raw).hexdigest()==meta['report_sha256'],'fixed accepted checkpoint report')
    report=json.loads(report_raw)
    need(report['passed'] is True and report['force_certificate_sha256']==SHA256
        and report['processed_nodes']==meta['processed_nodes']==25790
        and report['open_leaf_boxes']==meta['pending_leaves']==11
        and report['full_eight_point_exclusion_proved'] is False,'accepted incomplete checkpoint counts')
    names=set()
    for row in meta['arithmetic_inputs']:
        name=row['path'];p=PurePosixPath(name)
        need(not p.is_absolute() and '..' not in p.parts and ':' not in name and '\\' not in name
            and p.as_posix()==name and name not in names,'safe unique checkpoint arithmetic input')
        names.add(name)
        value=(ROOT/name).read_bytes()
        need(len(value)==row['bytes'] and hashlib.sha256(value).hexdigest()==row['sha256'],
            'unchanged accepted checkpoint arithmetic '+name)
    need(len(names)==meta['arithmetic_input_count']==50,'all checkpoint arithmetic inputs checked')
    return base,report,meta


def descendant(data,base):
    need(data.get('checkpoint_sha256')==SHA256,'explicit accepted checkpoint identity')
    for key in ['parent_sha256','attachment_sha256','precision_bits','scope','kernel_version','roots']:
        need(data[key]==base[key],'unchanged checkpoint model and root coverage '+key)
    nodes=data['nodes']
    need(len(nodes)>=len(base['nodes']),'checkpoint node prefix retained')
    for i,node in enumerate(base['nodes']):
        if node['kind']!='unknown':need(nodes[i]==node,'every accepted checkpoint proof step preserved')
