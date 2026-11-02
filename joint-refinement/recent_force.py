"""Identity and arithmetic binding for the already replayed 8be9f91 forest."""
from pathlib import Path,PurePosixPath
import gzip,hashlib,json

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
SHA256='173bf7f6ccf804a687afbb122af5b9130625f455ef2c104e1e3f2dbef390184f'
META_SHA256='ffb6037fd34a501005a7c63cdd41f8ee06863b76ac2e753a15524b467090d2f1'


def need(ok,label):
    if not ok:raise ArithmeticError(label)


def load_recent():
    raw=(HERE/'certificates/force_checkpoint_8be9f91.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest()==META_SHA256,'fixed accepted 8be9f91 metadata')
    meta=json.loads(raw)
    raw=(HERE/'certificates/force_checkpoint_8be9f91.json.gz').read_bytes()
    need(meta['certificate_sha256']==SHA256 and hashlib.sha256(raw).hexdigest()==SHA256,'fixed accepted 8be9f91 certificate')
    plain=gzip.decompress(raw)
    need(hashlib.sha256(plain).hexdigest()==meta['decompressed_sha256'],'accepted prefix decompressed identity')
    base=json.loads(plain)
    raw=(HERE/'certificates/force_checkpoint_8be9f91_report.json').read_bytes()
    need(hashlib.sha256(raw).hexdigest()==meta['report_sha256'],'fixed accepted 8be9f91 report')
    report=json.loads(raw)
    need(report['passed'] is True and report['force_certificate_sha256']==SHA256
        and report['processed_nodes']==meta['processed_nodes']==60791
        and report['open_leaf_boxes']==meta['pending_leaves']==30
        and report['full_eight_point_exclusion_proved'] is False,'accepted incomplete 8be9f91 prefix')
    names=set()
    for row in meta['arithmetic_inputs']:
        name=row['path'];p=PurePosixPath(name)
        need(not p.is_absolute() and '..' not in p.parts and ':' not in name and '\\' not in name
            and p.as_posix()==name and name not in names,'safe unique recent-prefix arithmetic input')
        names.add(name);raw=(ROOT/name).read_bytes()
        need(len(raw)==row['bytes'] and hashlib.sha256(raw).hexdigest()==row['sha256'],
            'unchanged accepted 8be9f91 arithmetic '+name)
    need(len(names)==meta['arithmetic_input_count']==61,'all recent-prefix arithmetic inputs checked')
    return base,report,meta


def descendant(data,base):
    need(data.get('accepted_prefix_sha256')==SHA256,'explicit accepted 8be9f91 prefix identity')
    for key in ['parent_sha256','attachment_sha256','precision_bits','scope','kernel_version','checkpoint_sha256','roots']:
        need(data[key]==base[key],'unchanged accepted-prefix model '+key)
    need(len(data['nodes'])>=len(base['nodes']),'complete accepted prefix retained')
    for i,node in enumerate(base['nodes']):
        if node['kind']!='unknown':need(data['nodes'][i]==node,'every accepted 8be9f91 step is unchanged')
