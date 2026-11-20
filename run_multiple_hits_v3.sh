#!/bin/bash
#SBATCH --job-name=REINVENT4_batch
#SBATCH --account=scripps-dept
#SBATCH --qos=scripps-dept
#SBATCH --partition=hpg-turin
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=reinvent_v3_%j.out
#SBATCH --error=reinvent_v3_%j.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
INPUT=${1:?Usage: sbatch run_multiple_hits_v3.sh input.txt [count] [output] [extra Python options...]}
N=${2:-1000}
OUT=${3:-REINVENT_V3_Results}
if (( $# >= 3 )); then shift 3; else shift "$#"; fi
/blue/mdisney/a.taghavi/reinvent4_env/bin/python - "$INPUT" "$N" "$OUT" "$@" <<'PY'
import csv, json, re, subprocess, sys
from pathlib import Path
from rdkit import Chem
source, target, root = Path(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3])
extra = sys.argv[4:]
if target < 1: raise SystemExit('Count must be positive')
if any(x.split('=')[0] in {'--out','--smiles','--input','--num-molecules','-n'} for x in extra):
    raise SystemExit('Extra options must not override input, output or count')
parents=[]
for line_no,line in enumerate(source.read_text().splitlines(),1):
    if not line.strip() or line.lstrip().startswith('#'): continue
    parts=line.split()
    if len(parts)==1: name=f'Hit_{len(parents)+1:04d}'; smi=parts[0]
    elif len(parts)==2: name,smi=parts
    else: raise SystemExit(f'Invalid input format at line {line_no}: expected SMILES or ID SMILES')
    safe=re.sub(r'[^A-Za-z0-9_.-]','_',name)
    folder=f'{len(parents)+1:04d}_{safe}' if len(parts)==2 else name
    parents.append((name,smi,folder))
if not parents: raise SystemExit('No inputs')
if root.exists() and any(root.iterdir()): raise SystemExit(f'Output directory not empty: {root}')
root.mkdir(parents=True,exist_ok=True)
base=['--model','priors/libinvent.prior','--device','cuda:0','--qed-min','0.7','--sa-max','4.0','--sc-max','3.5','--scscore-repo','./scscore','--scscore-model','scscore/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz','--max-rounds','30','--batch-size','256','--similarity-min','0.0']
(root/'batch_settings.json').write_text(json.dumps({'input':str(source),'target_per_input':target,'options':base+extra},indent=2))
all_rows=[]; summary=[]; unique={}; fields=[]
def write(path, rows, cols):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=cols); w.writeheader(); w.writerows(rows)
def merge():
    meta=['ParentID','ParentSMILES','ParentFolder','GeneratedID']
    write(root/'combined_generated.csv',all_rows,meta+fields)
    rows=[]
    for i,(smi, members) in enumerate(unique.items(),1):
        row=dict(members[0]); row['SMILES']=smi
        row['LibraryID']=f'LIB_{i:07d}'
        row['AllParentIDs']=json.dumps(list(dict.fromkeys(x['ParentID'] for x in members)))
        row['AllParentFolders']=json.dumps(list(dict.fromkeys(x['ParentFolder'] for x in members)))
        row['OccurrenceCount']=len(members); rows.append(row)
    write(root/'combined_unique_generated.csv',rows,['LibraryID','AllParentIDs','AllParentFolders','OccurrenceCount']+meta+fields)
    with (root/'combined_unique_generated.smi').open('w') as f:
        for row in rows: f.write(f"{row['SMILES']}\t{row['LibraryID']}\n")
    write(root/'overall_summary.csv',summary,['ParentID','ParentSMILES','ParentFolder','Requested','Accepted','Status','ReturnCode'])
failed=False
for name,smi,folder in parents:
    dest=root/folder
    print(f'Processing {name}: {smi}',flush=True)
    rc=subprocess.run([sys.executable,'scaffold_libinvent_v3.py','--smiles',smi,'--num-molecules',str(target),'--out',str(dest)]+base+extra).returncode
    count=0; path=dest/'generated.csv'
    if rc==0 and path.is_file():
        with path.open(newline='') as f:
            reader=csv.DictReader(f)
            for col in reader.fieldnames or []:
                if col not in fields: fields.append(col)
            for row in reader:
                mol=Chem.MolFromSmiles(row['SMILES'])
                if mol is None: raise ValueError(f'Invalid generated SMILES in {path}')
                key=Chem.MolToSmiles(mol,isomericSmiles=True)
                count+=1
                row.update(ParentID=name,ParentSMILES=smi,ParentFolder=folder,GeneratedID=f'{folder}_{count:06d}')
                all_rows.append(row); unique.setdefault(key,[]).append(row)
    status='failed' if rc or not path.is_file() else ('complete' if count>=target else 'shortfall')
    failed |= status=='failed'
    summary.append(dict(ParentID=name,ParentSMILES=smi,ParentFolder=folder,Requested=target,Accepted=count,Status=status,ReturnCode=rc))
    merge()
print(f'Combined: {len(all_rows)} records; {len(unique)} unique structures',flush=True)
if failed: sys.exit(1)
PY
