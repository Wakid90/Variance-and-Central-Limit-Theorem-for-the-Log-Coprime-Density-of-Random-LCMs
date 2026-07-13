.PHONY: reproduce paper test clean

reproduce:
	python code/reproduce_all.py

paper:
	cd paper && latexmk -pdf main.tex

test:
	pytest -q

clean:
	cd paper && latexmk -C
