# makefile for opt-syn

.PHONY: docs

docs:
	@$(MAKE) -C docs dirhtml
	 touch _build/dirhtml/.nojekyll
