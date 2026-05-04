PYTHON ?= python3

.PHONY: smoke validate

smoke:
	bash -n bin/ont_amplicon_vp1_by_folder.sh
	bash -n bin/vp1_pipeline_internal.sh
	$(PYTHON) -m py_compile bin/report_tables.py bin/report_html.py
	$(PYTHON) bin/report_tables.py --help >/dev/null
	$(PYTHON) bin/report_html.py --help >/dev/null

validate: smoke
	$(PYTHON) tests/validate_resources.py
