HTML_FILES := $(patsubst %.md, %.html, $(shell find . -name '*.md' \! -name 'navbar.md'))

all: $(HTML_FILES)

%.html: %.md navbar.html template.html
	pandoc '$<' -o '$@' --filter pandoc-include-code --template template.html --include-before navbar.html

navbar.html: navbar.md
	pandoc '$<' -o '$@'

.PHONY: clean
clean:
	$(RM) $(HTML_FILES) navbar.html

