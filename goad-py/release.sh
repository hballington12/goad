#!/bin/bash

# release.sh - Release management script for GOAD-PY
# Handles version bumping, tagging, and triggering CI releases

set -e

CURRENT_VERSION=$(grep '^version = ' pyproject.toml | cut -d'"' -f2)

echo "🚀 GOAD-PY Release Manager"
echo "=========================="
echo "Current version: $CURRENT_VERSION"
echo ""

# Parse arguments
ACTION=${1:-"help"}

show_help() {
    cat << EOF
Usage: $0 <command> [version]

Commands:
  patch         Bump patch version (0.2.0 -> 0.2.1)
  minor         Bump minor version (0.2.0 -> 0.3.0)  
  major         Bump major version (0.2.0 -> 1.0.0)
  custom <ver>  Set specific version (e.g., 0.2.1-rc1)
  tag           Create git tag for current version
  status        Show current version and git status
  help          Show this help

Examples:
  $0 patch      # Release 0.2.1
  $0 minor      # Release 0.3.0
  $0 custom 0.2.1-beta1  # Release beta
  $0 tag        # Tag current version

Release Process:
1. Run tests: ./test_wheels.sh
2. Bump version: $0 patch|minor|major
3. Create tag: $0 tag (triggers CI release)
4. Monitor: GitHub Actions will build and publish
EOF
}

bump_version() {
    local bump_type=$1
    local new_version
    
    IFS='.' read -ra VERSION_PARTS <<< "$CURRENT_VERSION"
    local major=${VERSION_PARTS[0]}
    local minor=${VERSION_PARTS[1]} 
    local patch=${VERSION_PARTS[2]}
    
    case $bump_type in
        "patch")
            new_version="$major.$minor.$((patch + 1))"
            ;;
        "minor")
            new_version="$major.$((minor + 1)).0"
            ;;
        "major")
            new_version="$((major + 1)).0.0"
            ;;
        "custom")
            new_version=$2
            if [ -z "$new_version" ]; then
                echo "❌ Error: Custom version required"
                echo "Usage: $0 custom 0.2.1-rc1"
                exit 1
            fi
            ;;
        *)
            echo "❌ Error: Invalid bump type: $bump_type"
            exit 1
            ;;
    esac
    
    echo "📝 Updating version: $CURRENT_VERSION -> $new_version"
    
    # Update pyproject.toml
    sed -i "s/version = \"$CURRENT_VERSION\"/version = \"$new_version\"/" pyproject.toml
    
    # Update Cargo.toml
    sed -i "s/version = \"$CURRENT_VERSION\"/version = \"$new_version\"/" Cargo.toml
    
    echo "✅ Version updated to $new_version"
    echo ""
    echo "📝 Next steps:"
    echo "1. Review changes: git diff"
    echo "2. Test build: ./build_and_test.sh"
    echo "3. Commit: git add . && git commit -m 'Bump version to $new_version'"
    echo "4. Tag: $0 tag"
}

create_tag() {
    local version=$(grep '^version = ' pyproject.toml | cut -d'"' -f2)
    local tag="v$version"
    
    # Check if working directory is clean
    if [ -n "$(git status --porcelain)" ]; then
        echo "⚠️  Warning: Working directory has uncommitted changes"
        echo "Commit your changes before tagging:"
        git status --short
        echo ""
        read -p "Continue anyway? (y/N): " confirm
        if [[ ! $confirm == [yY] ]]; then
            exit 1
        fi
    fi
    
    echo "🏷️  Creating tag: $tag"
    
    # Check if tag already exists
    if git tag -l | grep -q "^$tag$"; then
        echo "❌ Error: Tag $tag already exists"
        exit 1
    fi
    
    # Create annotated tag
    git tag -a "$tag" -m "Release $tag"
    
    echo "✅ Tag created: $tag"
    echo ""
    echo "🚀 To trigger release:"
    echo "git push origin $tag"
    echo ""
    echo "This will:"
    echo "- Trigger GitHub Actions CI"
    echo "- Build wheels for all platforms" 
    echo "- Run tests"
    echo "- Publish to PyPI automatically"
}

show_status() {
    echo "📊 Current Status:"
    echo "Version: $CURRENT_VERSION"
    echo "Git branch: $(git branch --show-current)"
    echo "Git status: $(git status --porcelain | wc -l) uncommitted changes"
    echo ""
    echo "Recent tags:"
    git tag -l | tail -5
}

case $ACTION in
    "patch"|"minor"|"major"|"custom")
        bump_version $ACTION $2
        ;;
    "tag")
        create_tag
        ;;
    "status")
        show_status
        ;;
    "help"|*)
        show_help
        ;;
esac