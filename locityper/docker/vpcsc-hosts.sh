# Sourced by non-interactive bash via BASH_ENV (dsub sidecar is bash -c).
# Docker overwrites /etc/hosts at container start.
_vpcsc_vip=199.36.153.4
for _vpcsc_h in \
    restricted.googleapis.com \
    storage.googleapis.com \
    oauth2.googleapis.com \
    www.googleapis.com \
    accounts.google.com \
    iamcredentials.googleapis.com
do
  if ! grep -q "[[:space:]]${_vpcsc_h}$" /etc/hosts 2>/dev/null; then
    echo "${_vpcsc_vip} ${_vpcsc_h}" >> /etc/hosts 2>/dev/null || true
  fi
done
unset _vpcsc_vip _vpcsc_h
