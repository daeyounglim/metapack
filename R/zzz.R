.onAttach<-function(libname, pkgname){
      # figure out year automatically (probably could be done more elegantly)
   date <- date()
   x <- regexpr("[0-9]{4}", date)
   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
   # echo output to screen
   packageStartupMessage("## Meta-Analysis (metapack)")
   packageStartupMessage("## Copyright (C) ", this.year,
      " Daeyoung Lim and Ming-Hui Chen", sep="")
}