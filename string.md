## AC自动机(拓扑优化)
ch:树边，转移边  
nex:回跳边
```
struct trie{
	int idx,ch[N][26],nex[N],sum[N],ind[N];
	vector<int>cnt[N];
	trie (){
		idx=0;
		memset(ind,0,sizeof(ind));
		memset(ch,0,sizeof(ch));
		memset(nex,0,sizeof(nex));
		memset(sum,0,sizeof(sum));
	}
	void insert(string s,int num){
		int p=0;
		for(int i=0;i<s.size();i++){
			int j=s[i]-'a';
			if(!ch[p][j])ch[p][j]=++idx;
			p=ch[p][j];
		}
		cnt[p].push_back(num);
	}
	void build(){
		queue<int>q;
		for(int i=0;i<26;i++){
			if(ch[0][i])q.push(ch[0][i]);
		}
		while(q.size()){
			int u=q.front();q.pop();
			for(int i=0;i<26;i++){
				int v=ch[u][i];
				if(v) nex[v]=ch[nex[u]][i],ind[nex[v]]++,q.push(v);
				else ch[u][i]=ch[nex[u]][i];
			}
		}
	}
	void query(string s,int n){
		for(int i=0,k=0;k<s.size();k++){
			i=ch[i][s[k]-'a'];
			sum[i]++;
		}
		vector<int>ans(n+1,0);
		queue<int>q;
		for(int i=1;i<=idx;i++){
			if(!ind[i])q.push(i);
		}
		while(q.size()){
			int u=q.front(),v=nex[u];q.pop();
			for(auto t:cnt[u]){
				ans[t]+=sum[u];
			}
			sum[v]+=sum[u];
			ind[v]--;
			if(ind[v]==0)q.push(v);
		}
		for(int i=1;i<=n;i++)cout<<ans[i]<<'\n';
	}
};
```
## 可持久化01trie
```
struct trie{
	int idx,cnt;
	int rt[N],ch[N*25][2],siz[N*25];
	void insert(int v){
		rt[++idx]=++cnt;
		int x=rt[idx-1];
		int y=rt[idx];
		for(int i=24;i>=0;i--){
			int j=v>>i&1;
			ch[y][!j]=ch[x][!j];//异位继承
			ch[y][j]=++cnt;//新位开点
			x=ch[x][j];
			y=ch[y][j];//走位
			siz[y]=siz[x]+1;//新位多1 
		}
	}
	int query(int x,int y,int v){//(x,y]左闭右开 
		int ans=0;
		x=rt[x],y=rt[y];
		for(int i=24;i>=0;i--){
			int j=v>>i&1;
			if(siz[ch[y][!j]]>siz[ch[x][!j]])
				x=ch[x][!j],y=ch[y][!j],ans|=1<<i;
			else{
				x=ch[x][j],y=ch[y][j];
			}
		}
		return ans;
	}
};
```
## KMP
```
void get_fail(string s,int fail[]){
	int j=0;
	fail[0]=0;
	for(int i=1;i<s.size();i++){
		while(j&&s[i]!=s[j])j=fail[j-1];
		if(s[i]==s[j])j++;
		fail[i]=j;
	}
}
```
```
string s,t;
cin>>s>>t;
get_fail(t,fail);
for(int i=0,j=0;i<s.size();i++){
	while(j&&s[i]!=t[j])j=fail[j-1];
	if(s[i]==t[j])j++;
	if(j==t.size()){
		cout<<i-t.size()+2<<'\n';
		j=fail[j-1];
	}
}
```
## Manacher
```
int C=0,R=0,p[N],ans;
string t,s="^#";
int n;
cin>>t;
n=t.size();
for(int i=0;i<n;i++){
	s+=t[i];
	s+="#";
}
s+="$";
for(int i=1;i<=2*n+1;i++){
	p[i]=i<R?min(R-i,p[2*C-i]):1;
	while(s[i+p[i]]==s[i-p[i]])p[i]++;
	if(i+p[i]>R)C=i,R=i+p[i];
	ans=max(ans,p[i]-1);//原回文串长度等于p[i]-1
}