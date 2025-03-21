import os
import subprocess
import shutil
from pathlib import Path

def upload_html_to_github(repo_url, analysis_output, local_repo_path, token, branch="gh-pages", commit_message="Deploy HTMLs to GitHub Pages"):
    try:
        cwd = os.getcwd()
        git_path = r"C:\\Users\\lwoods\\AppData\\Local\\GitHubDesktop\\app-3.4.9\\resources\\app\\git\\cmd\\git.exe"

        # Find the _data folder
        data_folder = Path(analysis_output) / (Path(analysis_output).name + "_data")
        if not data_folder.exists():
            print(f"Error: {data_folder} does not exist.")
            return

        print(f"Uploading folder: {data_folder}")

        # Ensure local repo exists
        if not Path(local_repo_path).exists():
            print(f"Cloning repository into {local_repo_path}...")
            subprocess.run([git_path, "clone", repo_url, local_repo_path], check=True)

        os.chdir(local_repo_path)

        # Set authenticated remote URL
        authenticated_repo_url = repo_url.replace("https://", f"https://{token}@")
        subprocess.run([git_path, "remote", "set-url", "origin", authenticated_repo_url], check=True)

        # Checkout the correct branch
        subprocess.run([git_path, "fetch", "origin"], check=True)
        subprocess.run([git_path, "checkout", branch], check=True)
        subprocess.run([git_path, "pull", "--rebase", "origin", branch], check=True)

        # Define destination for _data folder inside the repo
        dest_data_folder = Path(local_repo_path) / data_folder.name

        # Remove old _data folder in the repo
        if dest_data_folder.exists():
            shutil.rmtree(dest_data_folder)

        # Copy only .html and .png files, retaining structure
        exts = (".html", ".png")
        for root, _, files in os.walk(data_folder):
            for file in files:
                if file.endswith(exts):
                    src_file = Path(root) / file
                    rel_path = src_file.relative_to(data_folder)
                    dest_file = dest_data_folder / rel_path

                    # Ensure the directory exists
                    dest_file.parent.mkdir(parents=True, exist_ok=True)

                    # Copy the file
                    shutil.copy2(src_file, dest_file)
                    print(f"Copied: {src_file} -> {dest_file}")

        # Stage only the _data folder
        subprocess.run([git_path, "add", str(dest_data_folder)], check=True)

        # Check if there are staged changes before committing
        status_output = subprocess.run([git_path, "status", "--porcelain"], capture_output=True, text=True).stdout

        if status_output.strip():  # If there are changes
            subprocess.run([git_path, "commit", "-m", commit_message], check=True)
            subprocess.run([git_path, "push", "origin", branch], check=True)
            print("Changes pushed successfully.")
        else:
            print("No changes detected. Skipping commit & push.")


        # Push changes safely
        subprocess.run([git_path, "push", "origin", branch], check=True)

        print("HTML & PNG files successfully uploaded to GitHub Pages.")
        os.chdir(cwd)

    except subprocess.CalledProcessError as e:
        print(f"Git error: {e}")

    except Exception as e:
        print(f"Unexpected error: {e}")

