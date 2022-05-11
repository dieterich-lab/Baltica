pipeline {
    agent any
    stages {
        stage('Prepare venv') {
            steps {
                sh '''
                    python3 -m venv venv
                    . ./venv/bin/activate
                    pip install -r requirements.txt
                    pip install flake8 pytest pyfakefs
                '''
            }
          
        }
        stage('Linting') {
            steps {
                sh '''
                    . ./venv/bin/activate
                    flake8 baltica --max-line-length 140 --per-file-ignores baltica/parse_gffcompare_stats.py:F821
                '''
            }
        }
        stage('Testing') {
            steps {
                sh '''
                    . ./venv/bin/activate
                    export PYTHONPATH=$(pwd)
                    pytest --junitxml results.xml .tests
                '''
            }
        }
    }
    post {
        always {
            junit 'results.xml'
        }
    }
}
