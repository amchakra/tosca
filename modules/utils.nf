#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

def hiclipheader() {

    return """
    -----------------------------------------------------------------
    
    \033[0;34m           MMP""MM""YMM                                
               P'   MM   `7                                
                    MM  ,pW"Wq.  ,pP"Ybd  ,p6"bo   ,6"Yb.  
                    MM 6W'   `Wb 8I   `" 6M'  OO  8)   MM  
                    MM 8M     M8 `YMMMa. 8M        ,pm9MM  
                    MM YA.   ,A9 L.   I8 YM.    , 8M   MM  
                  .JMML.`Ybmd9'  M9mmmP'  YMbmd'  `Moo9^Yo.
    
    
    \033[0;36m${workflow.manifest.name} v${workflow.manifest.version}
    ${workflow.manifest.description}
    ${workflow.manifest.author}

    \033[0m-----------------------------------------------------------------
    """.stripIndent()

}
                                            
// MMP""MM""YMM                                
// P'   MM   `7                                
    //  MM  ,pW"Wq.  ,pP"Ybd  ,p6"bo   ,6"Yb.  
    //  MM 6W'   `Wb 8I   `" 6M'  OO  8)   MM  
    //  MM 8M     M8 `YMMMa. 8M        ,pm9MM  
    //  MM YA.   ,A9 L.   I8 YM.    , 8M   MM  
//    .JMML.`Ybmd9'  M9mmmP'  YMbmd'  `Moo9^Yo.
                                            
                                            
