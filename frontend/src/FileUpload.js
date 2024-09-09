import React, { useState } from 'react';
import axios from 'axios';

const FileUpload = () => {
    const [selectedFile, setSelectedFile] = useState(null);
    const [images, setImages] = useState({});

    const handleFileChange = (event) => {
        setSelectedFile(event.target.files[0]);
    };

    const handleUpload = async () => {
        if (!selectedFile) {
            alert('Please select a file.');
            return;
        }

        const formData = new FormData();
        formData.append('file', selectedFile);

        try {
            const response = await axios.post('http://localhost:5000/upload', formData, {
                headers: {
                    'Content-Type': 'multipart/form-data'
                }
            });

            setImages(response.data);

        } catch (error) {
            console.error('Error uploading file:', error);
        }
    };

    return (
        <div>
            <h1>mDNA Analysis Tool</h1>
            <input type="file" onChange={handleFileChange} />
            <button onClick={handleUpload}>Upload</button>

            <div>
                {Object.keys(images).length > 0 && (
                    <div>
                        <h2>Results</h2>
                        {images.gc_content && <img src={`http://localhost:5000/download/${images.gc_content.split("\\").pop()}`} alt="GC Content" />}
                        {images.sequence_length_distribution && <img src={`http://localhost:5000/download/${images.sequence_length_distribution.split("\\").pop()}`} alt="Sequence Length Distribution" />}
                        {images.codon_usage && <img src={`http://localhost:5000/download/${images.codon_usage.split("\\").pop()}`} alt="Codon Usage" />}
                        {images.conserved_regions && <img src={`http://localhost:5000/download/${images.conserved_regions.split("\\").pop()}`} alt="Conserved Regions" />}
                        {images.phylogenetic_tree && <img src={`http://localhost:5000/download/${images.phylogenetic_tree.split("\\").pop()}`} alt="Phylogenetic Tree" />}
                        
                        {images.gene_annotations && images.gene_annotations.map((annotation, index) => (
                            <img key={index} src={`http://localhost:5000/download/${annotation.split("\\").pop()}`} alt={`Gene Annotation ${index}`} />
                        ))}
                    </div>
                )}
            </div>
        </div>
    );
};

export default FileUpload;
