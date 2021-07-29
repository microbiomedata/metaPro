workflow myWorkflow {
    call myTask
}

task myTask {
    command {
        echo "hello world"
    }
    output {
        String out = read_string(stdout())
    }
}

task run_job_analysis {

    command {
        docker exec -it metapro_analysis-job_1 python3.8 ./metaProw/src/analysis_jobs/run_analysis_job.py
    }
    output {
        File out = stdout()
    }
    runtime {
        docker: 'metapro_analysis-job_1'
    }
}

workflow metaPro {
    call run_job_analysis
}