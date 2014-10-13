setClass("riboCoding", representation(CDS = "GRanges",
                                      hits = "array",
                                      unqHits = "array",
                                      replicates = "factor"))
setClass("riboData", representation(riboGR = "list",
                                    rnaGR = "list",
                                    replicates = "factor"))
