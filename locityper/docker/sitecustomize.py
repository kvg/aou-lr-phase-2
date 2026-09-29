# VPC-SC Batch VMs with --use-private-address often have no DNS for
# *.googleapis.com → restricted VIP. Docker regenerates /etc/hosts at start,
# so this runs on every python3 (dsub prepare/logging use python3 first).
_VIP = "199.36.153.4"
_HOSTS = (
    "restricted.googleapis.com",
    "storage.googleapis.com",
    "oauth2.googleapis.com",
    "www.googleapis.com",
    "accounts.google.com",
    "iamcredentials.googleapis.com",
)


def _ensure_restricted_vip_hosts() -> None:
    try:
        existing = open("/etc/hosts", encoding="utf-8").read()
    except OSError:
        return
    needed = [h for h in _HOSTS if h not in existing]
    if not needed:
        return
    try:
        with open("/etc/hosts", "a", encoding="utf-8") as fh:
            for host in needed:
                fh.write(f"{_VIP} {host}\n")
    except OSError:
        return


_ensure_restricted_vip_hosts()
