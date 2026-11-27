# ink

An implicit CAD kernel for agents, and the packages built on it.

| Folder | What |
|---|---|
| [kernel/](kernel/) | The engine, the `.ink` language and its standard library, and `ink`, the command line and MCP server. Start with [its README](kernel/README.md). |
| [sim/](sim/) | The simulation packages, built on the kernel as a third party's would be: studies, topology and dynamics, and an `ink` with them. [Its README](sim/README.md). |

Each is a Cargo workspace of its own: build and test it from its folder.

```sh
(cd kernel && cargo test)
(cd sim && cargo test)
```

[What changed, and upgrading from 0.1](kernel/CHANGELOG.md). [Contributing](kernel/CONTRIBUTING.md). MIT [license](LICENSE).
