import hashlib,json
from importlib import resources,metadata
import pytest

def test_installed_scientific_reference_hashes():
    root=resources.files('mumomac.resources')
    manifest=json.loads(root.joinpath('reference_version_manifest.json').read_text())
    for name,expected in {**manifest['resource_sha256'],**manifest['inherited_s16_resource_sha256']}.items():
        assert hashlib.sha256(root.joinpath(name).read_bytes()).hexdigest()==expected
    assert manifest['original_memberships']==229
    assert manifest['expanded_memberships']==820
    assert hashlib.sha256(root.joinpath('reference_registry_v1.json').read_bytes()).hexdigest()=='2c3ff81f43010caf11251053ec676d7d03ded3f941b9ac1660d836b5beeee515'

def test_installed_frozen_engines_and_versions():
    from mumomac_experimental.versions import validate_pinned_cores
    assert validate_pinned_cores()['integrity']=='PASS'
    assert metadata.version('mumomac-toolkit')=='0.3.0'
    assert metadata.version('mumomac-experimental')=='0.1.1rc1'

def test_no_private_developer_paths_in_distributed_manifests():
    for package,name in [('mumomac.resources','signatures_manifest.json'),('mumomac_experimental','resources/core_versions.json')]:
        text=resources.files(package).joinpath(name).read_text()
        assert '/home/' not in text
        assert '/mnt/' not in text

def test_missing_reference_resource_fails_clearly(tmp_path,monkeypatch):
    # An isolated subprocess in the CLI tests covers real package loading. Here
    # only the file location is replaced, retaining the original failure logic.
    import mumomac.reference as ref
    with monkeypatch.context() as m:
        m.setattr(ref.resources,'files',lambda name:tmp_path)
        with pytest.raises(FileNotFoundError,match='reference_version_manifest.json'):
            ref.load_reference()
