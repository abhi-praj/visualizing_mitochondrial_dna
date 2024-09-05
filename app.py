from flask import Flask, request, send_file, jsonify
import os
from werkzeug.utils import secure_filename
from GC_Analysis import display_gcc
# from ConservedRegions import conserved_regions_analysis
# from CodonUsageAnalysis import codon_usage_analysis
# from SequenceLengthDistribution import sequence_length_distribution
# from GeneAnnotation import gene_annotation_analysis
# from PhylogeneticTreeConstruction import phylogenetic_tree_construction

app = Flask(__name__)
UPLOAD_FOLDER = 'uploads'
RESULTS_FOLDER = 'results'

# Ensure upload and results directories exist
os.makedirs(UPLOAD_FOLDER, exist_ok=True)
os.makedirs(RESULTS_FOLDER, exist_ok=True)

# Allowed file types
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

        # Call all analysis functions here
        results = run_all_analyses(filepath)

        return jsonify(results)

    return jsonify({'error': 'File not allowed'}), 400


def run_all_analyses(filepath):
    """Run all DNA analyses and save results as image files."""
    results = {}

    # GC Content Analysis
    gc_image_path = os.path.join(RESULTS_FOLDER, 'gc_content.png')
    display_gcc(filepath, gc_image_path)
    results['gc_content'] = gc_image_path

    return results


@app.route('/download/<filename>')
def download_file(filename):
    """Allow users to download analysis result images."""
    return send_file(os.path.join(RESULTS_FOLDER, filename), as_attachment=True)


if __name__ == '__main__':
    app.run(debug=True)
