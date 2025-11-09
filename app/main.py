from fastapi import FastAPI, UploadFile, File, Form
from fastapi.responses import HTMLResponse
import io, base64

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from app.mutation_load import load_variant_file, compute_mutation_load

app = FastAPI()


@app.get("/", response_class=HTMLResponse)
def index():
    return """
    <html><body style="font-family:sans-serif;margin:40px">
      <h2>Mutation Load Visualizer</h2>
      <p>Upload 1–4 TSV/CSV files that contain at least chromosome, position, and alt-allele frequency.</p>
      <form action="/analyze" method="post" enctype="multipart/form-data">
        <p>File 1 (required): <input type="file" name="file1" required></p>
        <p>File 2 (optional): <input type="file" name="file2"></p>
        <p>File 3 (optional): <input type="file" name="file3"></p>
        <p>File 4 (optional): <input type="file" name="file4"></p>
        <p>Chromosome: <input type="text" name="chrom" value="V"></p>
        <p>Start: <input type="number" name="start" value="5292480"></p>
        <p>End: <input type="number" name="end" value="5293019"></p>
        <p><button type="submit">Plot</button></p>
      </form>
    </body></html>
    """


@app.post("/analyze", response_class=HTMLResponse)
async def analyze(
    file1: UploadFile = File(...),
    file2: UploadFile = File(None),
    file3: UploadFile = File(None),
    file4: UploadFile = File(None),
    chrom: str = Form(...),
    start: int = Form(...),
    end: int = Form(...),
):
    try:
        uploads = [file1, file2, file3, file4]
        dfs = []
        names = []

        for f in uploads:
            if f is None:
                continue
            content = await f.read()
            if not content:
                continue
            df = load_variant_file(content)
            dfs.append(df)
            names.append(f.filename or "uploaded.tsv")

        if not dfs:
            return HTMLResponse("<p style='color:red'>No valid files uploaded.</p>", status_code=400)

        fig, ax = plt.subplots(figsize=(8, 4))
        for i, df in enumerate(dfs):
            xs, ys = compute_mutation_load(df, chrom, start, end, bin_size=30)
            ax.plot(xs, ys, label=names[i], lw=2)

        ax.set_title(f"Mutation load {chrom}:{start}-{end}")
        ax.set_xlabel("Genomic position")
        ax.set_ylabel("Weighted score")
        ax.grid(alpha=0.4)
        ax.legend()

        buf = io.BytesIO()
        plt.tight_layout()
        fig.savefig(buf, format="png")
        plt.close(fig)
        img_b64 = base64.b64encode(buf.getvalue()).decode()

        return f"""
        <html><body style="font-family:sans-serif;margin:40px">
          <h2>Mutation-load plot</h2>
          <p>Files processed: {", ".join(names)}</p>
          <img src="data:image/png;base64,{img_b64}" />
          <p><a href="/">Back</a></p>
        </body></html>
        """

    except Exception as e:
        # THIS is what will show up in the browser instead of a plain 500
        return HTMLResponse(
            f"<h3 style='color:red'>Error while processing:</h3><pre>{e}</pre>",
            status_code=500,
        )
