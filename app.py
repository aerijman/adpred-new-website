from flask import Flask, request, render_template, url_for, send_from_directory
from bin.utils import *
from uuid import uuid4
from flask_mail import Mail, Message
import os, string

import resource, time

app = Flask(__name__)

app.config.update({
    'DEBUG': True,
    'MAIL_SERVER': 'smtp.gmail.com',
    'MAIL_PORT': 465,
    'MAIL_USE_TLS': False,
    'MAIL_USE_SSL': True,
    'MAIL_USERNAME': 'adprepredictor@gmail.com',
    'MAIL_PASSWORD': 'mejjEp-kekso2-cymqus',
})


mail = Mail(app)


@app.route('/', methods=['GET','POST'])
def index():

    if request.method == 'GET':
        return render_template('index.html', email_not_provided_msg='                   ')

    if request.method == 'POST':
        #print(request.form)

        #t1 = time.time()        

        email = request.form['email']
        if not re.match("[^@]+@[^@]+\.[^@]+",email):
            return render_template('index.html', email_not_provided_msg=" valid email must be provided!")
        
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

    plot = create_plot(adscore, sequence)

    #print(time.time()-t1, 'seconds')
    #print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)


    #print(sequence, struct, adscore, email)
    html_results = 'rendered_results/' + str(uuid4()) + '.html'


    # USE THIS TO SAVE THE RENDERED RESULTS AND SEND EMAIL WITH LINK WHEN IS READY. ALSO HAVE A CRON PROCESS THAT DELETES FILES OLDER THAN 24 HOURS.
    result_template = render_template("results.html", name='ADPred', plot=plot, csv_data=csv, csv_data_smooth=csv2, index='/index')
    with open(html_results,'w') as f:
        f.write(result_template)

    print('*********************',getBaseURL())
   

    ####################
    msg = Message(subject="ADpred results are ready", 
                  body="Your ADpred results are available at {}".format(getBaseURL() + '/' + html_results), 
                  sender="adprepredictor@gmail.com", 
                  recipients=[email])
    mail.send(msg)
    #####################

    return result_template
    #return render_template("results.html", name='ADPred', plot=plot, csv_data=csv, csv_data_smooth=csv2)

@app.route('/')
def Render(html):
    return render_template(html)

@app.route('/predictions/<filename>')
def download(filename):
    return send_from_directory('predictions', filename=filename) 


if __name__ == '__main__':
    app.run(port=5000, debug=True)
