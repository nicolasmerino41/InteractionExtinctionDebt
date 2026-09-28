"""Build the offline explorer after build_data.R. Uses only Python's standard library."""
from pathlib import Path
import json
import sys
folder = Path(__file__).resolve().parent
sys.path.insert(0, str(folder / 'renderer' / 'scripts'))
from render import export_html
data = json.loads((folder / 'datasets.json').read_text(encoding='utf-8'))
encoded = json.dumps(data, ensure_ascii=False, separators=(',', ':')).replace('</', '<\\/')
template = (folder / 'explorer-template.html').read_text(encoding='utf-8-sig')
fragment = folder / 'interaction-explorer.html'
fragment.write_text(template.replace('DATA_PLACEHOLDER', encoded), encoding='utf-8')
export_html(fragment, folder / 'index.html', title='Spatial interaction explorer', force=True)
print('Open', folder / 'index.html')
