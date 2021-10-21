include { checkm } from './process/checkm.nf'

workflow checkm_wf {
    take:   fastadir
    main:   
            checkm(fastadir)
}

