pipeline {
    stages {
        stage('prepare venv') {
            sh '''
                python3 -m venv venv
                . ./venv/bin/activate
                pip install -r requirements.txt
                pip install flake8 pytest pyfakefs
            '''
          
        }
        stage('linting') {
            sh '''
                . ./venv/bin/activate
                lake8 baltica --max-line-length 140 --per-file-ignores baltica/parse_gffcompare_stats.py:F821
            '''
        }
        stage('testing') {
            sh '''
                . ./venv/bin/activate
                export PYTHONPATH=$(pwd)
                pytest --junitxml results.xml .tests
            '''
        }
    }
    post {
        always {
            junit 'results.xml'
        }
    }
}
