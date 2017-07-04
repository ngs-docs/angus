# GitHub

From GitHub's about page: "**GitHub** is how people build software"

GitHub is a website meant to make the `git` version control system actually useable for humans.

The main page is: www.github.com

GitHub is great for:
 * Storing code
 * Contributing to other people's code
 * Getting help with other people's code
 * Keeping track of **issues**

## The Dashboard

Once you're logged in, GitHub will bring you to a page with several important sections.
 * On the left is your newsfeed, which tells you about updates to repositories or users you watch
 * In the top right is messages from GitHub
 * On the right are all your **repositories**

**Repositories** are locations where code (and other things) are stored and placed under **version control**. The repository is on one GitHub's servers.

**Version Control** is a system of keeping track of changes to a document.
It's like Track Changes in Microsoft Word, but better.
With version control, it is easy to "undo" work that may have been done weeks or even months ago.

The version control system GitHub uses is called `git`.

## Making Repositories

The main function of GitHub is to store repositories, so let's create one. To do so, look on the right side for the button that says "New repository". It is next to the section that lists your current repositories.

Click the button to create a new repository.

In the next window, enter a name into the box. You can also optionally enter a description in the box below. Additionally, you can initialize the repository with a `README.md` file. That means that when the repository is created, it will start by having a text file called `README.md`. The extension `md` stands for `Markdown` format that we saw earlier brifely in the [Jupyter Notebook lesson]() and we will cover more in depth in the [RMarkdown lesson]().

## Cloning Repositories

Congratulations! You now have a repository. To get it on your computer, you must now **clone** the repository. This will create a folder in your current directory and download all of the files -- as well as the entire history of the repository. It will also configure your computer to upload and download the repository.

## Commits

The way `git` implements version control is through something called a **commit**. A **commit** is a snapshot of a repository. It's similar to saving a file, except instead of one file we're making a commitment to the entire repository.

Making a commit is a two-step process:
 1. Add to the commit using `git add filename.txt`
 2. Do the commit using `git commit -m "helpful message"`

To upload the commit to GitHub:
 * `git push`
 * If that doesn't work, do `git push origin master`

To update your local copy *from* GitHub:
 * git pull

**WARNING: ALWAYS** `git pull` **BEFORE** you start working or git will complain.

## Issues

A useful feature of GitHub is the ability to keep track of issues.
Issues remind you of problems in your code.
You can also open issues on other people's code.
It's best to be as specific as possible and give as much information as possible to help the maintainers fix your problem.

You can cite issues in commit messages and they will show up in the issue. You do so by referencing the number of th issue. For example: `git commit -m "work on #1"`.
You can close issues through GitHub or through a commit message, by including the number in the message along with a word like "fix" or "fixed". For example: `git commit -m "fixed #1"`.

## Forking and Pull Requests

**Forking** a repository will copy it and create a repository in your account, under your control. It is completely separate from the repository it was forked from, though GitHub is aware of its source. To fork a repository, click the `Fork` button in the top right.

A **Pull Request** allows you to contribute the changes in your fork back to the original repository. Of course, this requires the permission of the owner of that repository.

## GitHub Education

GitHub understands the value of helping students learn their tool.
If you're a student, you can get free private repositories and many other professional services for free at https://education.github.com/








