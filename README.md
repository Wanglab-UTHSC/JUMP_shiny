# Running

The main repo is hosted on:

```
git clone git.genenetwork.org:/home/git/shared/source/jumpshiny
```

## Guix

On the server we can deploy the service as a Guix container. Basically install Guix and run `guix pull`:

```
guix pull # pull the latest version of guix
```

It takes a while. You only need to do this once!

Next, make sure the new Guix is in the path:

```
unset GUIX_PROFILE
~/.guix-profile/etc/profile
guix describe
```

Next, git clone guix-bioinformatics and run the script

```
git clone tux02.genenetwork.org:/home/git/public/guix-bioinformatics
source .guix-run
  (ignore warnings)
  Listening on http://127.0.0.1:3978
```

If you have a problem share the output on the console at any stage.

## nginx

For an example for proxy forwarding nginx, see [nginx.conf](etc/nginx.conf)
