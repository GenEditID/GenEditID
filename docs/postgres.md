##### Installing Postgres

Notes on installing Postgres on RHEL 7.

For this, I have installed Postgres 9.3 simply because the Clarity servers are running 9.3 and it makes Anne and my life easier to have the same version. Postgres is on 9.6 but I very much doubt we'll be using features in the later one that aren't in the earlier.

The shared instance is running on bioinf-srv003 until we get a dedicated VM.


#### Installing Packages

(As root)

    yum localinstall https://download.postgresql.org/pub/repos/yum/9.3/redhat/rhel-7-x86_64/pgdg-centos93-9.3-3.noarch.rpm
    yum install postgresql93-server postgresql93

    passwd postgres

(Set the password to something usable. It won't be used much.)

    su - postgres
    /usr/pgsql-9.3/bin/initdb
    exit

#### Configuring Access to the Database

(As postgres)

Edit `/var/lib/pgsql/9.3/data/pg_hba.conf` to add access to the (to be created) database `geneediting`. Add the following:

    # TYPE  DATABASE        USER        ADDRESS        METHOD
    local    geneediting    postgres                   ident
    local    geneediting    gene                       md5
    host     geneediting    gene        127.0.0.1/8    md5
    host     geneediting    gene        ::1/128        md5
    host     geneediting    gene        10.20.0.0/16   md5

Edit `/var/lib/pgsql/9.3/data/postgresql.conf` to change the `listen_addresses` to "`*`" (line 59).

(As root)

    systemctl enable postgresql-9.3
    systemctl start postgresql-9.3

(As postgres again)

    createdb geneediting
    createuser gene
    psql geneediting

(This is now in the Postgres client.)

    GRANT CREATE ON SCHEMA public TO "gene";
    GRANT USAGE ON SCHEMA public TO "gene";
    ALTER USER gene WITH PASSWORD 'gene';


#### Accessing the Database from Elsewhere

    psql -h bioinf-srv003.cri.camres.org -U gene geneediting

The YAML setting for accessing this database in the code is:

    DATABASE_URI: "postgresql://gene:gene@bioinf-srv003.cri.camres.org/geneediting"