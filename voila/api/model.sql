CREATE TABLE genome (
  row_id INTEGER PRIMARY KEY,
  name   VARCHAR NOT NULL
);


CREATE TABLE file_version (
  row_id INTEGER PRIMARY KEY,
  value  INTEGER NOT NULL
);


CREATE TABLE experiment (
  row_id INTEGER PRIMARY KEY,
  name   VARCHAR NOT NULL
);


CREATE TABLE alt_start (
  row_id     INTEGER PRIMARY KEY,
  gene_id    VARCHAR NOT NULL,
  coordinate INTEGER NOT NULL,
  FOREIGN KEY (gene_id) REFERENCES gene (row_id)
);


CREATE TABLE alt_end (
  row_id     INTEGER PRIMARY KEY,
  gene_id    VARCHAR NOT NULL,
  coordinate INTEGER NOT NULL,
  FOREIGN KEY (gene_id) REFERENCES gene (row_id)
);


CREATE TABLE junction_reads (
  row_id        INTEGER PRIMARY KEY,
  reads         INTEGER NOT NULL,
  experiment_id INTEGER NOT NULL,
  junction_id   INTEGER NOT NULL,
  FOREIGN KEY (junction_id) REFERENCES junction (row_id),
  FOREIGN KEY (experiment_id) REFERENCES experiment (row_id)
);


CREATE TABLE intron_retention_reads (
  row_id              INTEGER PRIMARY KEY,
  reads               INTEGER NOT NULL,
  experiment_id       VARCHAR NOT NULL,
  intron_retention_id INTEGER NOT NULL,
  FOREIGN KEY (intron_retention_id) REFERENCES intron_retention (row_id),
  FOREIGN KEY (experiment_id) REFERENCES experiment (row_id)
);


CREATE TABLE junction (
  row_id    INTEGER PRIMARY KEY,
  gene_id   VARCHAR NOT NULL,
  start     INTEGER NOT NULL,
  "end"     INTEGER NOT NULL,
  has_reads BOOLEAN,
  annotated BOOLEAN,
  FOREIGN KEY (gene_id) REFERENCES gene (row_id)
);


CREATE TABLE exon (
  row_id          INTEGER PRIMARY KEY,
  gene_id         VARCHAR NOT NULL,
  start           INTEGER NOT NULL,
  "end"           INTEGER NOT NULL,
  annotated_start INTEGER,
  annotated_end   INTEGER,
  annotated       BOOLEAN,
  FOREIGN KEY (gene_id) REFERENCES gene (row_id)
);


CREATE TABLE intron_retention (
  row_id    INTEGER PRIMARY KEY,
  gene_id   VARCHAR NOT NULL,
  start     INTEGER NOT NULL,
  "end"     INTEGER NOT NULL,
  has_reads BOOLEAN,
  annotated BOOLEAN,
  FOREIGN KEY (gene_id) REFERENCES gene (row_id)
);


CREATE TABLE gene (
  row_id     INTEGER PRIMARY KEY,
  id         VARCHAR NOT NULL,
  name       VARCHAR,
  strand     VARCHAR,
  chromosome VARCHAR
)

