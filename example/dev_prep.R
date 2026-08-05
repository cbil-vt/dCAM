# In dev mode, need to add Java path after devtools::load_all()

rJava::.jaddClassPath("./inst/java/")
rJava::.jaddClassPath("./inst/java/CornerDetect.jar")
# .jaddClassPath("./inst/java/lib/pj20150107.jar")
rJava::.jclassPath()

data <- matrix(c(0.1,0.2,1.0,0.0,0.0,0.5,0.3,
                 0.1,0.7,0.0,1.0,0.0,0.5,0.3,
                 0.8,0.1,0.0,0.0,1.0,0.0,0.4), nrow =3, byrow = TRUE)
topconv <- cornerSort(data, 3, 10)
