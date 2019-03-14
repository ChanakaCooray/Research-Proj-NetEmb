import git2json

fixture = open('temp/test2.txt')
commits = list(git2json.parse_commits(fixture.read()))
print(commits)