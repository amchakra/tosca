#!/usr/bin/env nextflow

// Specify DSL2
nextflow.enable.dsl=2

def hiclipheader() {

    return """
    -----------------------------------------------------------------
    
    \033[0;34m        __        _    ______  _____     _____  _______   
            [  |      (_) .' ___  ||_   _|   |_   _||_   __ \\  
            | |--.   __ / .'   \\_|  | |       | |    | |__) | 
            | .-. | [  || |         | |   _   | |    |  ___/  
            | | | |  | |\\ `.___.'\\ _| |__/ | _| |_  _| |_     
            [___]|__][___]`.____ .'|________||_____||_____|  
    
    
    \033[0;36m${workflow.manifest.name} v${workflow.manifest.version}
    ${workflow.manifest.description}
    ${workflow.manifest.author}

    \033[0m-----------------------------------------------------------------
    """.stripIndent()

}