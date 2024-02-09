import tempfile

import pyarrow as pa


def loose(record_batches, schema, tmpdir=None):
    """Loose an arrow :) Stream compressed text files into temporary arrow files

    polars + arrow = very fast reading and processing
    """
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp()

    arrowpath = tempfile.NamedTemporaryFile(dir=tmpdir, delete=False)

    with pa.OSFile(arrowpath.name, "wb") as sink:
        with pa.RecordBatchFileWriter(sink=sink, schema=schema) as writer:
            for batch in record_batches:
                writer.write(batch)

    return arrowpath
