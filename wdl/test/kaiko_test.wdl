version 1.0

import "../workflow/kaiko.wdl" as kaiko

workflow kaiko_run_test {
    input  {
        File raw_file
        File kaiko_config
        File mgf_file
    }

    call kaiko.run {
        input:
            raw_file = raw_file,
            kaiko_config = kaiko_config,
            mgf_file = mgf_file
    }
}
