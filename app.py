from flask import Flask, request, render_template, url_for, send_from_directory
from bin.utils import *
from uuid import uuid4

app = Flask(__name__)

@app.route('/', methods=['GET','POST'])
def index():

    if request.method == 'GET':
        return render_template('index.html', email_not_provided_msg='                   ')

    if request.method == 'POST':
        #print(request.form)

        email = request.form['email']
        if email == '':
            return render_template('index.html', email_not_provided_msg=" email must be provided!")
        
        # check email exists
        #if mail.does-not-exists()...:
        #   return render_template('error.html', errorType=error_no_mail)

        # User provides only protein ID
        if request.form['Sequence'] == '':
            protId = request.form['protIdName']
            sequence = identifier2fasta(protId)
            
            if sequence == -1:
                return render_template('error.html', errorType=error_uniprot)

        # User provides sequence
        else:
            sequence = request.form['Sequence']        

    # remove fasta header
    sequence = clean_input(sequence)
    
    # run psipred
    randName = str(uuid4())
    with open(randName, 'w') as f:
        f.write(sequence)

    struct = get_psipred(randName, email)
    adscore, csv, csv2 = predict_full(sequence, struct, randName)

    print(sequence, struct, adscore, email)

    plot = create_plot(adscore, sequence)
        
    return render_template("results.html", name='ADPred', plot=plot, csv_data=csv, csv_data_smooth=csv2)
    

@app.route('/predictions/<filename>')
def download(filename):
    return send_from_directory('predictions', filename=filename) 


if __name__ == '__main__':
    app.run(port=5000, debug=True)
