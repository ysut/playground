MAKEFLAGS += --no-builtin-rules
MAKEFLAGS += --warn-undefined-variables
SHELL := /bin/bash

.PHONY: up
.DEFAULT_GOAL := help

help: 
	@echo "Usage: make [target]"
	@echo ""
	@echo "Targets:"
	@echo "  up           Start the application"
	@echo "  exec         Connect to the mariadb container"
	@echo "  destruction  !! Destroy the application !!"
	@echo "  help         Show this help message"

up:
	docker compose up --build --detach --wait \
	&& docker compose exec mariadb mariadb \
		--defaults-extra-file=/home/sql.cnf -e "SHOW DATABASES;"

restore:
	docker compose exec mariadb mariadb \
		--defaults-extra-file=/home/sql.cnf hgmd < /home/hgmd.sql

exec:
	docker compose exec mariadb bash

destruction:
	docker compose down --rmi all --volumes --remove-orphans
	sudo rm -rf mariadb_data/*