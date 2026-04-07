.PHONY: bootstrap docker-build

bootstrap:
	bash bin/bootstrap.sh

docker-build:
	docker build -t mr-gwas-pipeline:latest .
