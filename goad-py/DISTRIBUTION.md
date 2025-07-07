# GOAD-PY Distribution Guide

## 🚀 Quick Release Process

```bash
# 1. Test everything works
./build_and_test.sh
./test_wheels.sh

# 2. Bump version and create release
./release.sh patch          # or minor/major
git add . && git commit -m "Bump version to X.Y.Z"
./release.sh tag
git push origin vX.Y.Z      # Triggers automatic release

# 3. Monitor release
# GitHub Actions will automatically:
# - Build for all platforms
# - Run tests
# - Publish to PyPI
```

## 📋 Pre-Release Checklist

### ✅ Local Testing
- [ ] `./build_and_test.sh` passes
- [ ] `./test_wheels.sh` passes
- [ ] Examples work (`simple_example.py`, `multiproblem_example.py`)
- [ ] Type stubs are up to date (`goad_py.pyi`)

### 🧪 TestPyPI Testing (Optional)
```bash
# Manual TestPyPI upload
./publish_test.sh

# Test installation from TestPyPI
pip install --index-url https://test.pypi.org/simple/ goad-py
```

### 📝 Release Preparation
- [ ] Update version in `pyproject.toml` and `Cargo.toml`
- [ ] Update `README.md` if needed
- [ ] Commit all changes
- [ ] Create git tag with `./release.sh tag`

## 🏗️ Platform Coverage

Our CI builds wheels for:

### Linux (manylinux)
- x86_64, x86, aarch64, armv7, s390x, ppc64le

### Linux (musllinux) 
- x86_64, x86, aarch64, armv7

### Windows
- x64, x86

### macOS
- x86_64 (Intel), aarch64 (Apple Silicon)

### Testing
- Ubuntu, Windows, macOS
- Python 3.8, 3.9, 3.10, 3.11, 3.12

## 🔐 Security Features

- **Trusted Publishing**: No API tokens needed (uses OIDC)
- **Build Attestation**: Cryptographic proof of build provenance
- **Artifact Signing**: All releases are signed
- **Supply Chain Security**: Full CI/CD audit trail

## 📚 Distribution Workflow

### Automatic (Recommended)
1. **Push tag** → GitHub Actions triggers
2. **Build wheels** for all platforms
3. **Run tests** on all platforms
4. **Generate attestations** for security
5. **Publish to PyPI** automatically

### Manual (Development)
1. **TestPyPI**: `./publish_test.sh`
2. **Local build**: `./build_and_test.sh`
3. **Wheel testing**: `./test_wheels.sh`

## 🛠️ Scripts Reference

| Script | Purpose |
|--------|---------|
| `build_and_test.sh` | Local development build |
| `build_wheels_local.sh` | Build with cibuildwheel |
| `test_wheels.sh` | Validate built wheels |
| `publish_test.sh` | Upload to TestPyPI |
| `release.sh` | Version management & tagging |

## 🔧 Setup Requirements

### For Contributors
- Rust toolchain
- Python 3.8+
- Git

### For CI/CD (One-time setup)
1. **GitHub Repository**: Enable Actions
2. **PyPI Account**: Create trusted publisher
   - Project: `goad-py`
   - Owner: `hballington12` 
   - Repository: `goad`
   - Workflow: `CI.yml`
   - Environment: `pypi`
3. **TestPyPI** (optional): Same as above for testing

### Trusted Publishing Setup
1. Go to PyPI → Account Settings → Publishing
2. Add GitHub publisher:
   - Repository: `hballington12/goad`
   - Workflow: `CI.yml`
   - Environment: `pypi`

## 📊 Release Monitoring

### GitHub Actions
- **View builds**: `Actions` tab in GitHub
- **Download artifacts**: Available for 90 days
- **View logs**: Detailed build information

### PyPI
- **Package page**: https://pypi.org/project/goad-py/
- **Download stats**: Available on PyPI
- **Version history**: All releases tracked

## 🐛 Troubleshooting

### Build Failures
- Check Rust compilation errors
- Verify all platforms build
- Review test failures

### Upload Failures  
- Check trusted publishing setup
- Verify tag format (`vX.Y.Z`)
- Review PyPI project settings

### Version Conflicts
- Ensure version bump before release
- Check for existing tags
- Verify semver format

## 📈 Success Metrics

A successful release should:
- ✅ Build for 15+ platform combinations
- ✅ Pass tests on 3 operating systems  
- ✅ Generate signed attestations
- ✅ Upload to PyPI within 30 minutes
- ✅ Be installable via `pip install goad-py`