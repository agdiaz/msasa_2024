run-all-experiments-1a:
	nextflow run main.nf --input_dir ./datasets/BB11004 -with-docker
	nextflow run main.nf --input_dir ./datasets/BB11013 -with-docker
	# nextflow run main.nf --input_dir ./datasets/BB30003 -with-docker
run-all-experiments-1b:
	nextflow run main.nf --input_dir ./datasets/BB40010 -with-docker
	nextflow run main.nf --input_dir ./datasets/BB40014 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBA0011 -with-docker

run-all-experiments-2a:
	nextflow run main.nf --input_dir ./datasets/BBA0030 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBA0065 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBA0117 -with-docker
run-all-experiments-2b:
	nextflow run main.nf --input_dir ./datasets/BBA0142 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBA0192 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBS12014 -with-docker

run-all-experiments-3a:
	nextflow run main.nf --input_dir ./datasets/BBS20002 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBS20007 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBS20009 -with-docker
run-all-experiments-3b:
	nextflow run main.nf --input_dir ./datasets/BBS30001 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBS30003 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBS30015 -with-docker

run-all-experiments-4a:
	nextflow run main.nf --input_dir ./datasets/BBS30016 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBS30017 -with-docker
	nextflow run main.nf --input_dir ./datasets/BBS50005 -with-docker
run-all-experiments-4b:
	nextflow run main.nf --input_dir ./datasets/BOX122 -with-docker
	nextflow run main.nf --input_dir ./datasets/BOX132 -with-docker
	nextflow run main.nf --input_dir ./datasets/BOX212 -with-docker

run-all:
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BB11004 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BB11013 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BB30003 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BB40010 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BB40014 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBA0011 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBA0030 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBA0065 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBA0117 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBA0142 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBA0192 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBS12014 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBS20002 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBS20007 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBS20009 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBS30001 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBS30003 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBS30015 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBS30016 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBS30017 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BBS50005 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BOX122 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BOX132 -with-docker
	NFX_OPTS="-Xms=1g -Xmx=8g" NXF_ANSI_LOG=false nextflow run main.nf --input_dir ./datasets/BOX212 -with-docker

test-run:
	NFX_OPTS="-Xms=1g -Xmx=8g" nextflow run main.nf --input_dir ./datasets/BB11013 -with-docker -with-conda -with-report reports/BB11013.html -with-timeline timeline/BB11013.html


run-experiment:
	# 1. Coincidences
	echo "Coincidences Strategy" > $(INPUT_FASTA_FILE).msasa.timing.log
	date +"%Y-%m-%d %H:%M:%S" >> $(INPUT_FASTA_FILE).msasa.timing.log

	/usr/bin/time -p python msasa_cli.py $(INPUT_FASTA_FILE) $(INPUT_FASTA_FILE).msasa.1.coincidences.aln \
		--log_file $(INPUT_FASTA_FILE).msasa.1.coincidences.log \
		--extend \
		--temp 10 \
		--cooling_rate 0.99 \
		--min_temp 0.001 \
		--max_no_changes 500 \
		--quality_function coincidences 2>> $(INPUT_FASTA_FILE).msasa.timing.log

	# 2. Identity
	echo "Identity Strategy" >> $(INPUT_FASTA_FILE).msasa.timing.log
	date +"%Y-%m-%d %H:%M:%S" >> $(INPUT_FASTA_FILE).msasa.timing.log

	/usr/bin/time -p python msasa_cli.py $(INPUT_FASTA_FILE) $(INPUT_FASTA_FILE).msasa.2.identity.aln \
		--log_file $(INPUT_FASTA_FILE).msasa.2.identity.log \
		--extend \
		--temp 5 \
		--cooling_rate 0.99 \
		--min_temp 0.0001 \
		--max_no_changes 500 \
		--quality_function identity 2>> $(INPUT_FASTA_FILE).msasa.timing.log

	# 3. Similarity
	echo "Similarity Strategy (Blosum62)" >> $(INPUT_FASTA_FILE).msasa.timing.log
	date +"%Y-%m-%d %H:%M:%S" >> $(INPUT_FASTA_FILE).msasa.timing.log

	/usr/bin/time -p python msasa_cli.py $(INPUT_FASTA_FILE) $(INPUT_FASTA_FILE).msasa.3.similarity.aln \
		--log_file $(INPUT_FASTA_FILE).msasa.3.similarity.log \
		--extend \
		--temp 5 \
		--cooling_rate 0.99 \
		--min_temp 0.01 \
		--max_no_changes 500 \
		--quality_function similarity 2>> $(INPUT_FASTA_FILE).msasa.timing.log

	# 4. Global strategy
	echo "Global Strategy" >> $(INPUT_FASTA_FILE).msasa.timing.log
	date +"%Y-%m-%d %H:%M:%S" >> $(INPUT_FASTA_FILE).msasa.timing.log

	/usr/bin/time -p python msasa_cli.py $(INPUT_FASTA_FILE) $(INPUT_FASTA_FILE).msasa.4.global.aln \
		--log_file $(INPUT_FASTA_FILE).msasa.4.global.log \
		--extend \
		--temp 10 \
		--cooling_rate 0.99 \
		--min_temp 0.01 \
		--max_no_changes 500 \
		--quality_function global 2>> $(INPUT_FASTA_FILE).msasa.timing.log

	# 5. Local Strategy
	echo "Local Strategy" >> $(INPUT_FASTA_FILE).msasa.timing.log
	date +"%Y-%m-%d %H:%M:%S" >> $(INPUT_FASTA_FILE).msasa.timing.log

	/usr/bin/time -p python msasa_cli.py $(INPUT_FASTA_FILE) $(INPUT_FASTA_FILE).msasa.5.local.aln \
		--log_file $(INPUT_FASTA_FILE).msasa.5.local.log \
		--extend \
		--temp 10 \
		--cooling_rate 0.99 \
		--min_temp 0.01 \
		--max_no_changes 500 \
		--quality_function local 2>> $(INPUT_FASTA_FILE).msasa.timing.log
