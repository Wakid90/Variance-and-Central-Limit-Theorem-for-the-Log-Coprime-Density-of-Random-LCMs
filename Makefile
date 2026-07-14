.PHONY: all repro paper test clean

all: repro paper

## Regenerate all figures, tables, and machine-readable outputs
repro:
	python code/reproduce_all.py

## Compile the manuscript
paper:
	cd paper && latexmk -pdf main.tex

## Run the numerical regression tests
test:
	pytest -q

clean:
	cd paper && latexmk -C
	rm -rf .pytest_cache **/__pycache__
