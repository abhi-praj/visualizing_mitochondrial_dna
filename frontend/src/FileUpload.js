import React, { useState } from 'react';
import axios from 'axios';

const FileUpload = () => {
    const [selectedFile, setSelectedFile] = useState(null);
    const [results, setResults] = useState({});

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

            setResults(response.data);

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
                {Object.keys(results).length > 0 && (
                    <div>
                        <h2>Results</h2>

                        {/* GC Content HTML */}
                        {results.gc_content && (
                            <div>
                                <h3>GC Content</h3>
                                <iframe
                                    src={`http://localhost:5000/download/${results.gc_content.split("\\").pop()}`}
                                    title="GC Content"
                                    width="100%"
                                    height="500px"
                                ></iframe>
                            </div>
                        )}

                        {/* Sequence Length Distribution HTML */}
                        {results.sequence_length_distribution && (
                            <div>
                                <h3>Sequence Length Distribution</h3>
                                <iframe
                                    src={`http://localhost:5000/download/${results.sequence_length_distribution.split("\\").pop()}`}
                                    title="Sequence Length Distribution"
                                    width="100%"
                                    height="500px"
                                ></iframe>
                            </div>
                        )}

                        {/* Codon Usage HTML */}
                        {results.codon_usage && (
                            <div>
                                <h3>Codon Usage</h3>
                                <iframe
                                    src={`http://localhost:5000/download/${results.codon_usage.split("\\").pop()}`}
                                    title="Codon Usage"
                                    width="100%"
                                    height="500px"
                                ></iframe>
                            </div>
                        )}

                        {/* Conserved Regions HTML */}
                        {results.conserved_regions && (
                            <div>
                                <h3>Conserved Regions</h3>
                                <iframe
                                    src={`http://localhost:5000/download/${results.conserved_regions.split("\\").pop()}`}
                                    title="Conserved Regions"
                                    width="100%"
                                    height="500px"
                                ></iframe>
                            </div>
                        )}
                    </div>
                )}
            </div>
        </div>
    );
};

export default FileUpload;
