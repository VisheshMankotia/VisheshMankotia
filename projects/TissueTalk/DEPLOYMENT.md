# Deploy TissueTalk with Streamlit Community Cloud

This project currently lives in the `projects/TissueTalk` subfolder of the `VisheshMankotia` repository.

1. Push the latest changes to the repository's `main` branch.
2. Sign in to [Streamlit Community Cloud](https://share.streamlit.io/) with GitHub.
3. Choose **Create app** and select repository `VisheshMankotia/VisheshMankotia`.
4. Select branch `main`.
5. Set the app file path to `projects/TissueTalk/app.py`.
6. In Advanced settings, select Python 3.11 (or another version compatible with the dependencies).
7. Deploy.

Streamlit should install dependencies from the app's `requirements.txt` in the same directory. If the platform does not detect the nested requirements file automatically, configure the deployment to use `projects/TissueTalk/requirements.txt` where supported, or place a root-level deployment requirements file.

## Before sharing publicly
- Test the app with the simulated demo dataset.
- Test with a small, non-sensitive AnnData file.
- Confirm the app's resource use with the intended dataset size.
- Do not upload confidential or restricted research data.
- Add a deployed URL to the project README only after deployment succeeds.
