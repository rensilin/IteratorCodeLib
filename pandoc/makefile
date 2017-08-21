html:index.html
	@echo index.html

index.html:../iteratorCodeLib.md
	pandoc --mathjax --template html.latex -s -N -o "$@" "$<"

kkke.me:index.html
	scp "$<" kkke.me:html

pdf:iteratorCodeLibrary.pdf
	@echo iteratorCodeLibrary.pdf

iteratorCodeLibrary.pdf:../iteratorCodeLib.md
	pandoc --latex-engine=xelatex --template pdf.tex -s --toc -N -o "$@" "$<"

