# Mitochondrial DNA Analysis & Visualization Platform

[![Python](https://img.shields.io/badge/Python-3.8+-blue.svg)](https://python.org)
[![React](https://img.shields.io/badge/React-18.3.1-blue.svg)](https://reactjs.org)
[![Flask](https://img.shields.io/badge/Flask-2.0+-green.svg)](https://flask.palletsprojects.com)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

A web-based platform for comprehensive mitochondrial DNA (mtDNA) analysis and visualization.

## Key Features

### Advanced DNA Analysis
- **GC Content Analysis**: Comprehensive analysis of guanine-cytosine content distribution across mitochondrial genomes
- **Sequence Length Distribution**: Statistical analysis and visualization of sequence length variations
- **Codon Usage Analysis**: Detailed examination of codon usage patterns and amino acid frequencies
- **Conserved Regions Detection**: Identification and visualization of evolutionarily conserved genomic regions

### Interactive Visualizations
- **Plotly-Powered Charts**: High-quality, interactive visualizations that can be zoomed, panned, and exported
- **Real-time Processing**: Fast analysis of large mitochondrial DNA datasets
- **Professional Reports**: Publication-ready HTML reports with embedded interactive charts

## Installation & Setup

### Prerequisites
- Python 3.8 or higher
- Node.js 14 or higher
- npm or yarn package manager

### How to Run

For the fastest setup, use our automated quick start script:

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/visualizing_mitochondrial_dna.git
   cd visualizing_mitochondrial_dna
   ```

2. **Run the quick start script**
   ```bash
   python start_project.py
   ```
   
   This script will:
   - Check system requirements
   - Install all dependencies automatically
   - Start both backend and frontend servers
   - Open your browser to the application

### Manual Setup

#### Backend Setup

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/visualizing_mitochondrial_dna.git
   cd visualizing_mitochondrial_dna
   ```

2. **Install Python dependencies**
   ```bash
   cd backend
   pip install -r requirements.txt
   ```

3. **Start the Flask server**
   ```bash
   python app.py
   ```
   The backend will be available at `http://localhost:5000`

#### Frontend Setup

1. **Navigate to frontend directory**
   ```bash
   cd frontend
   ```

2. **Install dependencies**
   ```bash
   npm install
   ```

3. **Start the React development server**
   ```bash
   npm start
   ```
   The frontend will be available at `http://localhost:3000`

## Usage Guide

### Basic Workflow

1. **Upload Data**: Select a FASTA file containing mitochondrial DNA sequences
2. **Automatic Analysis**: The system automatically runs four comprehensive analyses:
   - GC content distribution analysis
   - Sequence length distribution profiling
   - Codon usage pattern analysis
   - Conserved regions identification
3. **View Results**: Interactive visualizations are generated and displayed in real-time
4. **Export Reports**: Download publication-ready HTML reports with embedded charts

### Supported File Formats
- **FASTA (.fasta)**: Standard format for nucleotide sequences
- **Multi-FASTA**: Multiple sequences in a single file for comparative analysis

### Sample Data
The repository includes sample mitochondrial genome data from the National Library of Medicine for testing and demonstration purposes.

## Analysis Capabilities

### GC Content Analysis
- Calculates guanine-cytosine content percentage for each sequence
- Generates distribution histograms showing GC content variation
- Identifies sequences with unusual base composition patterns

### Sequence Length Distribution
- Analyzes sequence length variations across the dataset
- Creates statistical summaries of genome size distributions
- Helps identify incomplete or truncated sequences

### Codon Usage Analysis
- Translates DNA sequences to amino acid sequences
- Calculates codon usage frequencies and patterns
- Generates amino acid frequency distributions
- Useful for studying translational efficiency and evolutionary pressure

### Conserved Regions Detection
- Performs multiple sequence alignment
- Calculates conservation scores for each genomic position
- Visualizes conservation patterns across evolutionary time
- Identifies functionally important genomic regions

---
