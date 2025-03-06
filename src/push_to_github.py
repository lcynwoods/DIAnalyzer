"""Send to GitHub"""

import os
import subprocess
from pathlib import Path

import os
import subprocess
from pathlib import Path

def upload_html_to_github(repo_url, folder_path, token, branch="gh-pages", commit_message="Deploy HTMLs to GitHub Pages"):
    try:
        cwd = os.getcwd()
        # Update repo_url to include the token for authentication
        authenticated_repo_url = repo_url.replace("https://", f"https://{token}@")

        os.chdir(folder_path)
        git_path = r"C:\\Users\\lwoods\\AppData\\Local\\GitHubDesktop\\app-3.4.9\\resources\\app\\git\\cmd\\git.exe"

        # Initialize a new git repository (if it's not already a git repository)
        if not (Path(folder_path) / ".git").exists():
            subprocess.run([git_path, "init"], check=True)

        # Check if the remote 'origin' already exists
        try:
            existing_remotes = subprocess.check_output([git_path, "remote", "-v"], text=True).strip()
            if "origin" in existing_remotes:
                # Update the remote URL if it already exists
                subprocess.run([git_path, "remote", "set-url", "origin", authenticated_repo_url], check=True)
            else:
                # Add the remote repository
                subprocess.run([git_path, "remote", "add", "origin", authenticated_repo_url], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Failed to check or set the remote: {e}")

        # Checkout the branch for GitHub Pages
        print(1)
        subprocess.run([git_path, "checkout", "-B", branch], check=True)

        # Collect all .html files with a path length limit
        max_path_length = 260
        html_files = [
            os.path.relpath(os.path.join(root, file), folder_path)
            for root, _, files in os.walk(folder_path)
            for file in files if file.endswith(".html") and len(os.path.join(root, file)) < max_path_length
        ]

        ext = [".html", ".png"]

        viz_files = [
            os.path.relpath(os.path.join(root, file), folder_path)
            for root, _, files in os.walk(folder_path)
            for file in files if file.endswith(tuple(ext)) and len(os.path.join(root, file)) < max_path_length
        ]

        if not html_files:
            print("No HTML files found or all file paths exceed the maximum length. Exiting without making changes.")
            return

        print(3)
        # Remove previously staged files to avoid uploading non-HTML files
                # Check if there are any staged files before attempting to remove them
        try:
            staged_files = subprocess.check_output([git_path, "ls-files"], text=True).strip()
            if staged_files:
                # Remove previously staged files to avoid uploading non-HTML files
                subprocess.run([git_path, "rm", "-r", "--cached", "."], check=True)
            else:
                print("No staged files to remove.")
        except subprocess.CalledProcessError as e:
            print(f"Error checking staged files: {e}")

        print(4)
        # Add only HTML files to staging
        for relative_path in viz_files:
            subprocess.run([git_path, "add", relative_path], check=True)

        # Commit the changes
        print(5)
        subprocess.run([git_path, "commit", "-m", commit_message], check=True)

        # Push to the GitHub Pages branch
        subprocess.run([git_path, "push", "-u", "origin", branch, "--force"], check=True)

        print("HTML files successfully uploaded to GitHub Pages.")
        os.chdir(cwd)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred during a git operation: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")