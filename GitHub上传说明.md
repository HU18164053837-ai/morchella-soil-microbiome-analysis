# GitHub 一键上传说明

1. 先在 GitHub 网站新建一个**空仓库**，不要勾选自动添加 README、LICENSE 或 `.gitignore`。
2. 推荐仓库名：`morchella-soil-microbiome-analysis`。
3. 如果论文尚未正式发表，建议先设为 **Private**；投稿或数据公开时再改为 Public。
4. 双击仓库根目录中的 `一键检查并上传到GitHub.bat`。
5. 输入空仓库地址，例如 `https://github.com/你的用户名/morchella-soil-microbiome-analysis.git`。
6. 第一次推送时，GitHub 可能要求浏览器登录或 Personal Access Token。

脚本会在推送前检查大文件、原始序列、文稿、密钥、个人账号路径和硬编码本机路径。检查失败时不会推送。
