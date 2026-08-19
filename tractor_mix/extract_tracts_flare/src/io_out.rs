use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

use flate2::write::GzEncoder;
use flate2::Compression;

const BUF: usize = 1 << 20;

pub enum OutWriter {
    Plain(BufWriter<File>),
    Gzip(BufWriter<GzEncoder<BufWriter<File>>>),
}

impl OutWriter {
    pub fn create(path: &Path, gzip: bool) -> io::Result<Self> {
        let file = File::create(path)?;
        let inner = BufWriter::with_capacity(BUF, file);
        if gzip {
            let gz = GzEncoder::new(inner, Compression::new(6));
            Ok(Self::Gzip(BufWriter::with_capacity(BUF, gz)))
        } else {
            Ok(Self::Plain(inner))
        }
    }

    pub fn write_all(&mut self, buf: &[u8]) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.write_all(buf),
            Self::Gzip(w) => w.write_all(buf),
        }
    }

    pub fn finish(self) -> io::Result<()> {
        match self {
            Self::Plain(mut w) => w.flush(),
            Self::Gzip(w) => {
                let gz = w.into_inner()?;
                gz.finish()?.flush()?;
                Ok(())
            }
        }
    }
}

impl Write for OutWriter {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Self::Plain(w) => w.write(buf),
            Self::Gzip(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Self::Plain(w) => w.flush(),
            Self::Gzip(w) => w.flush(),
        }
    }
}
