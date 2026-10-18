"""Integrity adapter for the single exact frozen-core snapshot."""
import json,hashlib
from importlib import resources,metadata
def validate_pinned_cores():
 manifest=json.loads(resources.files('mumomac_experimental').joinpath('resources/core_versions.json').read_text())
 for row in manifest['files']:
  parts=row['packaged_path'].split('/')[1:];root=resources.files(parts[0]);data=root.joinpath(*parts[1:]).read_bytes()
  if hashlib.sha256(data).hexdigest()!=row['sha256']:raise RuntimeError('Pinned core differs: '+row['packaged_path'])
 if metadata.version('mumomac-toolkit')!='0.3.0':raise RuntimeError('Frozen formal Toolkit0.3.0 required.')
 return dict(versions=manifest['versions'],n_frozen_files=len(manifest['files']),integrity='PASS')
