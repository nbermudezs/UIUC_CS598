
## Development

### Pre-requisites

The following dependencies are assumed to be installed and working:
- Python 3
- virtualenv

### Setup
After cloning this repository please follow the steps below:
- run `virtualenv env --python=$(which python3)` to create a new virtual environment
- activate the newly created env `source env/bin/activate`
- install dependencies by running `pip install -r requirements.txt`
- download the data necessary `./utils/download_data.sh`
- preprocess the data `python data/preprocess.py`