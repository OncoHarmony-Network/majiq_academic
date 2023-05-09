CREATE TABLE genome (
  id   INTEGER NOT NULL,
  name VARCHAR,
  PRIMARY KEY (id)
);

CREATE TABLE file_version (
  id    INTEGER NOT NULL,
  value INTEGER,
  PRIMARY KEY (id)
);

CREATE TABLE experiment (
  name VARCHAR NOT NULL,
  PRIMARY KEY (name)
);

CREATE TABLE alt_start (
  gene_id    VARCHAR NOT NULL,
  coordinate INTEGER NOT NULL,
  PRIMARY KEY (gene_id, coordinate),
  FOREIGN KEY (gene_id) REFERENCES gene (id)
);

CREATE TABLE alt_end (
  gene_id    VARCHAR NOT NULL,
  coordinate INTEGER NOT NULL,
  PRIMARY KEY (gene_id, coordinate),
  FOREIGN KEY (gene_id) REFERENCES gene (id)
);

CREATE TABLE junction_reads (
  reads            INTEGER NOT NULL,
  experiment_name  VARCHAR NOT NULL,
  junction_gene_id VARCHAR NOT NULL,
  junction_start   INTEGER NOT NULL,
  junction_end     INTEGER NOT NULL,
  PRIMARY KEY (experiment_name, junction_gene_id, junction_start, junction_end),
  FOREIGN KEY (junction_gene_id, junction_start, junction_end) REFERENCES junction (gene_id, start, "end"),
  FOREIGN KEY (experiment_name) REFERENCES experiment (name)
);

CREATE TABLE intron_retention_reads (
  reads                    INTEGER NOT NULL,
  experiment_name          VARCHAR NOT NULL,
  intron_retention_gene_id VARCHAR NOT NULL,
  intron_retention_start   INTEGER NOT NULL,
  intron_retention_end     INTEGER NOT NULL,
  PRIMARY KEY (experiment_name, intron_retention_gene_id, intron_retention_start, intron_retention_end),
  FOREIGN KEY (intron_retention_gene_id, intron_retention_start, intron_retention_end) REFERENCES intron_retention (gene_id, start, "end"),
  FOREIGN KEY (experiment_name) REFERENCES experiment (name)
);

CREATE TABLE junction (
  gene_id   VARCHAR NOT NULL,
  start     INTEGER NOT NULL,
  "end"     INTEGER NOT NULL,
  has_reads BOOLEAN,
  annotated BOOLEAN,
  is_simplified BOOLEAN,
  is_constitutive BOOLEAN,
  has_flag BOOLEAN,
  PRIMARY KEY (gene_id, start, "end"),
  FOREIGN KEY (gene_id) REFERENCES gene (id),
  CHECK (has_reads IN (0, 1)),
  CHECK (annotated IN (0, 1)),
  CHECK (is_constitutive IN (0, 1)),
  CHECK (has_flag IN (0, 1))
);

CREATE TABLE exon (
  gene_id         VARCHAR NOT NULL,
  start           INTEGER NOT NULL,
  "end"           INTEGER NOT NULL,
  annotated_start INTEGER,
  annotated_end   INTEGER,
  annotated       BOOLEAN,
  PRIMARY KEY (gene_id, start, "end"),
  FOREIGN KEY (gene_id) REFERENCES gene (id),
  CHECK (annotated IN (0, 1))
);

CREATE TABLE intron_retention (
  gene_id   VARCHAR NOT NULL,
  start     INTEGER NOT NULL,
  "end"     INTEGER NOT NULL,
  has_reads BOOLEAN,
  annotated BOOLEAN,
  is_simplified BOOLEAN,
  is_constitutive BOOLEAN,
  has_flag BOOLEAN,
  PRIMARY KEY (gene_id, start, "end"),
  FOREIGN KEY (gene_id) REFERENCES gene (id),
  CHECK (has_reads IN (0, 1)),
  CHECK (annotated IN (0, 1)),
  CHECK (is_constitutive IN (0, 1)),
  CHECK (has_flag IN (0, 1))
);

CREATE TABLE  gene_overlap (
  gene_id_1 VARCHAR NOT NULL,
  gene_id_2 VARCHAR NOT NULL,
  PRIMARY KEY (gene_id_1, gene_id_2)
);

CREATE TRIGGER remove_dup BEFORE INSERT
ON gene_overlap
BEGIN
   DELETE FROM gene_overlap WHERE gene_id_1 = new.gene_id_1 and gene_id_2 = new.gene_id_2;
   DELETE FROM gene_overlap WHERE gene_id_1 = new.gene_id_2 and gene_id_2 = new.gene_id_1;
END;


CREATE TABLE gene (
  id         VARCHAR NOT NULL,
  name       VARCHAR,
  strand     VARCHAR,
  chromosome VARCHAR,
  PRIMARY KEY (id)
);

CREATE INDEX `alt_end_p` ON `alt_end` ( `gene_id`, `coordinate` );
CREATE INDEX `alt_start_p` ON `alt_start` ( `gene_id`, `coordinate` );
CREATE INDEX `exon_p` ON `exon` ( `gene_id`, `start`, `end` );
CREATE INDEX `experiment_p` ON `experiment` ( `name` );
CREATE INDEX `file_version_p` ON `file_version` ( `id` );
CREATE INDEX `gene_overlap_p` ON `gene_overlap` ( `gene_id_1`, `gene_id_2` );
CREATE INDEX `gene_p` ON `gene` ( `id` );
CREATE INDEX `genome_p` ON `genome` ( `id` );
CREATE INDEX `intron_retention_p` ON `intron_retention` ( `gene_id`, `start`, `end` );
CREATE INDEX `intron_retention_reads_p` ON `intron_retention_reads` ( `intron_retention_gene_id`, `intron_retention_start`, `intron_retention_end` , `experiment_name`);
CREATE INDEX `intron_retention_reads_p2` ON `intron_retention_reads` (`intron_retention_gene_id`);
CREATE INDEX `junction_p` ON `junction` ( `gene_id`, `start`, `end` );
CREATE INDEX `junction_reads_p` ON `junction_reads` ( `junction_gene_id`, `junction_start`, `junction_end`, `experiment_name`);
CREATE INDEX `junction_reads_p2` ON `junction_reads` ( `junction_gene_id` );

