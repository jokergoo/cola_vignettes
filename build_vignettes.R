library(knitr)
library(markdown)

markdown2html = function(rmd) {
	md = gsub("rmd$", "md", rmd, ignore.case = TRUE)
	html = gsub("rmd$", "html", rmd, ignore.case = TRUE)
	knit(rmd, md)
	markdownToHTML(md, html)
}

markdown2html("cola_quick.Rmd")
markdown2html("cola.Rmd")
markdown2html("functional_enrichment.Rmd")
markdown2html("predict.Rmd")
markdown2html("work_with_big_datasets.Rmd")
markdown2html("compare_partitions.Rmd")
markdown2html("hierarchical.Rmd")
