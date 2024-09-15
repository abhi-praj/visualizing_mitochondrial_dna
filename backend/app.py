from flask import Flask, request, send_file, jsonify
import os
from flask_cors import CORS
from werkzeug.utils import secure_filename
from GC_Analysis import display_gcc
from SequenceLengthDistribution import sequence_length_distribution
from CodonUsageAnalysis import analyze_codon_usage
from ConservedRegions import conserved_regions

app = Flask(__name__)
CORS(app, resources={r"/*": {"origins": "*"}})
UPLOAD_FOLDER = 'uploads'
RESULTS_FOLDER = 'results'

os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)

ALLOWED_EXTENSIONS = {'fasta'}

def allowed_file(filename):
    """Check if uploaded file has an allowed extension."""
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in ALLOWED_EXTENSIONS

@app.route('/')
def index():
    """Home route to welcome the user."""
    return "Welcome to the DNA Analysis Tool"

@app.route('/upload', methods=['POST'])
def upload_file():
    """Handle file upload and run analyses."""
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    if file and allowed_file(file.filename):
        filename = secure_filename(file.filename)
        filepath = os.path.join(UPLOAD_FOLDER, filename)
        file.save(filepath)

        results = run_all_analyses(filepath)

        return jsonify(results)

    return jsonify({'error': 'File not allowed'}), 400

def run_all_analyses(filepath):
    """Run all DNA analyses and save results as HTML files."""
    results = {}

    gc_html_path = os.path.join(RESULTS_FOLDER, 'gc_content.html')
    display_gcc(filepath, gc_html_path)
    results['gc_content'] = gc_html_path

    seq_length_html_path = os.path.join(RESULTS_FOLDER, 'sequence_length_distribution.html')
    sequence_length_distribution(filepath, seq_length_html_path)
    results['sequence_length_distribution'] = seq_length_html_path

    codon_usage_html_path = os.path.join(RESULTS_FOLDER, 'codon_usage.html')
    analyze_codon_usage(filepath, codon_usage_html_path)
    results['codon_usage'] = codon_usage_html_path

    conserved_regions_html_path = os.path.join(RESULTS_FOLDER, 'conserved_regions.html')
    conserved_regions(filepath, conserved_regions_html_path)
    results['conserved_regions'] = conserved_regions_html_path

    return results

@app.route('/download/<filename>')
def download_file(filename):
    """Allow users to download analysis result files."""
    file_path = os.path.join(RESULTS_FOLDER, filename)
    
    if filename.endswith('.html'):
        return send_file(file_path, as_attachment=False, mimetype='text/html')
    
    return send_file(file_path, as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)
