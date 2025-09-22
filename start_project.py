#!/usr/bin/env python3
"""
Quick start script for the Mitochondrial DNA Analysis Platform
This script helps users set up and run both backend and frontend services.
"""

import subprocess
import sys
import os
import time
import webbrowser
from pathlib import Path

def check_python_version():
    """Check if Python version is compatible."""
    if sys.version_info < (3, 8):
        print("Error: Python 3.8 or higher is required.")
        print(f"Current version: {sys.version}")
        return False
    print(f"Python version: {sys.version.split()[0]}")
    return True

def check_node_version():
    """Check if Node.js is installed."""
    try:
        result = subprocess.run(['node', '--version'], capture_output=True, text=True)
        if result.returncode == 0:
            print(f"Node.js version: {result.stdout.strip()}")
            return True
    except FileNotFoundError:
        pass
    print("❌ Error: Node.js is not installed or not in PATH.")
    print("Please install Node.js from https://nodejs.org/")
    return False

def install_backend_dependencies():
    """Install Python dependencies for the backend."""
    print("\n🔧 Installing backend dependencies...")
    backend_dir = Path("backend")
    requirements_file = backend_dir / "requirements.txt"
    
    if not requirements_file.exists():
        print("requirements.txt not found in backend directory")
        return False
    
    try:
        subprocess.run([sys.executable, "-m", "pip", "install", "-r", str(requirements_file)], 
                      check=True, cwd=backend_dir)
        print("Backend dependencies installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to install backend dependencies: {e}")
        return False

def install_frontend_dependencies():
    """Install Node.js dependencies for the frontend."""
    print("\nInstalling frontend dependencies...")
    frontend_dir = Path("frontend")
    package_json = frontend_dir / "package.json"
    
    if not package_json.exists():
        print("package.json not found in frontend directory")
        return False
    
    try:
        subprocess.run(["npm", "install"], check=True, cwd=frontend_dir)
        print("Frontend dependencies installed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"Failed to install frontend dependencies: {e}")
        return False

def start_backend():
    """Start the Flask backend server."""
    print("\nStarting backend server...")
    backend_dir = Path("backend")
    
    try:
        # Start backend in a separate process
        backend_process = subprocess.Popen([sys.executable, "app.py"], 
                                         cwd=backend_dir)
        print("Backend server started on http://localhost:5000")
        return backend_process
    except Exception as e:
        print(f"Failed to start backend server: {e}")
        return None

def start_frontend():
    """Start the React frontend server."""
    print("\nStarting frontend server...")
    frontend_dir = Path("frontend")
    
    try:
        # Start frontend in a separate process
        frontend_process = subprocess.Popen(["npm", "start"], 
                                          cwd=frontend_dir)
        print("Frontend server started on http://localhost:3000")
        return frontend_process
    except Exception as e:
        print(f"Failed to start frontend server: {e}")
        return None

def main():
    """Main function to set up and start the application."""
    print("🧬 Mitochondrial DNA Analysis Platform - Quick Start")
    print("=" * 55)
    
    # Check system requirements
    if not check_python_version():
        return
    
    if not check_node_version():
        return
    
    # Install dependencies
    if not install_backend_dependencies():
        return
    
    if not install_frontend_dependencies():
        return
    
    # Start servers
    print("\nStarting servers...")
    backend_process = start_backend()
    if not backend_process:
        return
    
    # Wait a moment for backend to start
    time.sleep(3)
    
    frontend_process = start_frontend()
    if not frontend_process:
        backend_process.terminate()
        return
    
    print("\n" + "=" * 55)
    print(" Application started successfully!")
    print("\n Access the application:")
    print("   Frontend: http://localhost:3000")
    print("   Backend API: http://localhost:5000")
    print("\n Tips:")
    print("   - Upload a FASTA file to start analyzing mitochondrial DNA")
    print("   - Check the 'sample_genome.fasta' file for example data")
    print("   - Press Ctrl+C to stop both servers")
    
    # Open browser after a short delay
    time.sleep(5)
    try:
        webbrowser.open("http://localhost:3000")
        print("Browser opened automatically")
    except:
        print("Please open http://localhost:3000 in your browser")
    
    try:
        # Keep the script running
        frontend_process.wait()
    except KeyboardInterrupt:
        print("\n\nShutting down servers...")
        frontend_process.terminate()
        backend_process.terminate()
        print("Servers stopped successfully")

if __name__ == "__main__":
    main()
